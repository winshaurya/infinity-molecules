import json
import uuid
import os
import sys
from supabase import create_client, Client
from dotenv import load_dotenv
load_dotenv()

sys.path.append(os.path.join(os.path.dirname(__file__), 'backend'))
from core_logic import generate_functionalized_isomers

SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")
supabase: Client = create_client(SUPABASE_URL, SUPABASE_KEY)

def get_user_id(request):
    auth_header = request.headers.get('Authorization')
    if not auth_header or not auth_header.startswith('Bearer '):
        return None
    token = auth_header.split(' ')[1]
    try:
        user_response = supabase.auth.get_user(token)
        return user_response.user.id
    except:
        return None

def handler(request):
    try:
        user_id = get_user_id(request)
        if not user_id:
            return {'error': 'Unauthorized'}, 401

        method = request.method

        if method == 'POST':
            data = request.get_json()
            job_id = data.get('job_id')
            molecules_count = data.get('molecules_count', 0)
            download_format = data.get('download_format', 'csv')

            return download_molecules(user_id, job_id, molecules_count, download_format)

        return {'error': 'Method not allowed'}, 405

    except Exception as e:
        return {'error': str(e)}, 500

def download_molecules(user_id, job_id, molecules_count, download_format):
    # Get job
    response = supabase.table('jobs').select('*').eq('id', job_id).eq('user_id', user_id).execute()
    if not response.data:
        return {'error': 'Job not found'}, 404

    job = response.data[0]
    if job['status'] != 'completed':
        return {'error': 'Job not completed'}, 400

    # Allow any molecules_count - we'll clamp it below

    # Get user
    user_response = supabase.table('users').select('*').eq('id', user_id).execute()
    if not user_response.data:
        return {'error': 'User not found'}, 404

    user = user_response.data[0]

    # Check tier
    tier = 'free'
    if user['is_fullaccess']:
        tier = 'fullaccess'
    elif user['credits'] >= 15:
        tier = 'paid'

    if download_format == 'all' and tier != 'fullaccess':
        return {'error': 'Full batch download requires fullaccess tier'}, 403

    # Clamp molecules_count to available molecules
    actual_molecules_count = min(molecules_count, job['total_molecules'])

    # Calculate credits based on actual molecules to download
    credits_required = (actual_molecules_count + 999) // 1000  # Ceiling division

    if user['credits'] < credits_required:
        return {'error': f'Insufficient credits. Required: {credits_required}, Available: {user["credits"]}'}, 402

    # Deduct credits
    new_credits = user['credits'] - credits_required
    supabase.table('users').update({'credits': new_credits}).eq('id', user_id).execute()

    # Add credit history
    supabase.table('credit_history').insert({
        'user_id': user_id,
        'amount': -credits_required,
        'reason': f'Download {actual_molecules_count} molecules',
        'description': f'Downloaded {actual_molecules_count} molecules in {download_format} format'
    }).execute()

    # Regenerate molecules (since not stored)
    params = json.loads(job['parameters'])
    smiles_list = generate_functionalized_isomers(
        n_carbons=params["carbon_count"],
        functional_groups=params["functional_groups"],
        n_double_bonds=params["double_bonds"],
        n_triple_bonds=params["triple_bonds"],
        n_rings=params["rings"],
        carbon_types=params["carbon_types"]
    )

    # Take actual count (clamped to available)
    selected_smiles = smiles_list[:actual_molecules_count]

    # Return SMILES data to frontend for client-side processing
    data = json.dumps({
        'smiles': selected_smiles,
        'format': download_format,
        'job_id': job_id
    })
    content_type = 'application/json'
    filename = f'molecules_{job_id}.json'

    return {
        'success': True,
        'credits_used': credits_required,
        'remaining_credits': new_credits,
        'actual_molecules_count': actual_molecules_count,
        'download_id': str(uuid.uuid4()),
        'data': data,
        'content_type': content_type,
        'filename': filename,
        'message': f'Download prepared. {credits_required} credits deducted for {actual_molecules_count} molecules.'
    }
