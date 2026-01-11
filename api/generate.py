import json
import uuid
import time
import datetime
import os
import sys
import signal
import gc
import threading
from functools import partial

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from supabase import create_client, Client
from dotenv import load_dotenv
load_dotenv()

# Import core logic
sys.path.append(os.path.join(os.path.dirname(__file__), 'backend'))
from core_logic import generate_functionalized_isomers, validate_structure_possibility

# Supabase setup
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")
supabase: Client = create_client(SUPABASE_URL, SUPABASE_KEY)

def handler(request):
    try:
        # Get token from Authorization header
        auth_header = request.headers.get('Authorization')
        if not auth_header or not auth_header.startswith('Bearer '):
            return {'error': 'Unauthorized'}, 401

        token = auth_header.split(' ')[1]

        # Verify token and get user
        try:
            user_response = supabase.auth.get_user(token)
            user_id = user_response.user.id
        except Exception as e:
            return {'error': 'Invalid token'}, 401

        # Parse request body
        data = request.get_json()
        carbon_count = data.get('carbon_count', 0)
        double_bonds = data.get('double_bonds', 0)
        triple_bonds = data.get('triple_bonds', 0)
        rings = data.get('rings', 0)
        carbon_types = data.get('carbon_types', ["primary", "secondary", "tertiary"])
        functional_groups = data.get('functional_groups', [])

        # Validate
        is_valid, error = validate_structure_possibility(
            carbon_count, functional_groups, double_bonds, triple_bonds, carbon_types, rings
        )
        if not is_valid:
            return {'error': error}, 400

        # Create job in DB with processing status
        job_id = str(uuid.uuid4())
        job_data = {
            'id': job_id,
            'user_id': user_id,
            'parameters': json.dumps(data),
            'status': 'processing',
            'total_molecules': 0  # Will be updated after generation
        }
        supabase.table('jobs').insert(job_data).execute()

        # Generate molecules (this is the time-consuming part)
        smiles_list = generate_functionalized_isomers(
            n_carbons=carbon_count,
            functional_groups=functional_groups,
            n_double_bonds=double_bonds,
            n_triple_bonds=triple_bonds,
            n_rings=rings,
            carbon_types=carbon_types
        )

        # Update job with completed status and actual molecule count
        supabase.table('jobs').update({
            'status': 'completed',
            'total_molecules': len(smiles_list)
        }).eq('id', job_id).execute()

        return {
            'job_id': job_id,
            'status': 'completed',
            'total_molecules': len(smiles_list),
            'message': 'Generation completed'
        }

    except Exception as e:
        return {'error': str(e)}, 500
