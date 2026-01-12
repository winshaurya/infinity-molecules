from fastapi import FastAPI, HTTPException, Depends, Request, Header
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field
from typing import List, Dict, Any
from fastapi.responses import HTMLResponse
import uuid
import time
from datetime import datetime
import os
import json
import sys
import os
import io
import zipfile
import base64
import gc
import asyncio
import functools
from concurrent.futures import ProcessPoolExecutor
from supabase import create_client, Client
from dotenv import load_dotenv
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from core_logic import generate_functionalized_isomers, validate_structure_possibility, get_functional_group_type, get_element_counts, rdMolDescriptors
import traceback

# lazy-import core logic & RDKit with graceful handling so the process doesn't crash
STARTUP_ERROR = None
try:
    from core_logic import generate_functionalized_isomers, validate_structure_possibility, get_functional_group_type, get_element_counts, rdMolDescriptors
    from rdkit import Chem
except Exception:
    STARTUP_ERROR = traceback.format_exc()
    # Export placeholders so module imports succeed; endpoints will return clear error
    generate_functionalized_isomers = None
    validate_structure_possibility = None
    get_functional_group_type = None
    get_element_counts = None
    rdMolDescriptors = None
    Chem = None

# Load environment variables
load_dotenv()

app = FastAPI(title="Chemistry SaaS API")

# CORS Configuration
FRONTEND_URL = os.getenv("FRONTEND_URL", "http://localhost:5173")
BACKEND_URL = os.getenv("BACKEND_URL", "http://localhost:8000")

app.add_middleware(
    CORSMiddleware,
    allow_origins=[FRONTEND_URL],  # Allow frontend domain
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Generation runtime safeguards
GENERATION_TIMEOUT_SECONDS = int(os.getenv("GENERATION_TIMEOUT_SECONDS", "180"))
GENERATION_MAX_WORKERS = int(os.getenv("GENERATION_MAX_WORKERS", "1"))

PROCESS_POOL = ProcessPoolExecutor(max_workers=GENERATION_MAX_WORKERS)

# Supabase setup
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_SERVICE_ROLE_KEY = os.getenv("SUPABASE_SERVICE_ROLE_KEY")
SUPABASE_ANON_KEY = os.getenv("SUPABASE_ANON_KEY")

print(f"Loading SUPABASE_URL: {SUPABASE_URL}")
print(f"Loading SUPABASE_SERVICE_ROLE_KEY: {SUPABASE_SERVICE_ROLE_KEY[:20] if SUPABASE_SERVICE_ROLE_KEY else 'Not set'}...")
print(f"Loading SUPABASE_ANON_KEY: {SUPABASE_ANON_KEY[:20] if SUPABASE_ANON_KEY else 'Not set'}...")

# Use service role key for server-side operations (bypasses RLS)
supabase_key = SUPABASE_SERVICE_ROLE_KEY

if not SUPABASE_URL or not supabase_key:
    raise ValueError("SUPABASE_URL and either SUPABASE_SERVICE_ROLE_KEY or SUPABASE_ANON_KEY must be set")

try:
    supabase: Client = create_client(SUPABASE_URL, supabase_key)
    print("✅ Supabase client initialized successfully")
except Exception as e:
    print(f"❌ Supabase initialization failed: {e}")
    raise

# Models
class JobRequest(BaseModel):
    carbon_count: int = Field(ge=0, le=20)
    double_bonds: int = Field(ge=0, le=10)
    triple_bonds: int = Field(ge=0, le=10)
    rings: int = Field(ge=0, le=10)
    carbon_types: List[str] = Field(default=["primary", "secondary", "tertiary"])
    functional_groups: List[str] = Field(max_length=50)

class DownloadRequest(BaseModel):
    job_id: str
    molecules_count: int = Field(gt=0)
    download_format: str = Field(pattern="^(csv|molsdf)$")

class CreditRefillRequest(BaseModel):
    amount: int = Field(gt=0, description="Amount of credits to add")
    description: str = Field(default="Manual credit refill")

# Database helper functions
def get_user_by_id(user_id: str):
    """Get user data by ID"""
    try:
        response = supabase.table('users').select('*').eq('id', user_id).execute()
        if response.data:
            user = response.data[0]
            return {
                'id': user['id'],
                'email': user['email'],
                'credits': user['credits'],
                'is_fullaccess': user['is_fullaccess'],
                'created_at': user['created_at'],
                'updated_at': user['updated_at']
            }
        return None
    except Exception as e:
        print(f"Error getting user {user_id}: {e}")
        raise

def update_user_credits(user_id: str, new_credits: int):
    """Update user credits"""
    try:
        supabase.table('users').update({
            'credits': new_credits,
            'updated_at': 'now()'
        }).eq('id', user_id).execute()
    except Exception as e:
        print(f"Error updating credits for user {user_id}: {e}")
        raise

def create_job(job_data: dict):
    """Create a new job"""
    try:
        supabase.table('jobs').insert(job_data).execute()
    except Exception as e:
        print(f"Error creating job: {e}")
        raise

def update_job_status(job_id: str, status: str, total_molecules: int = None, completed_at: str = None):
    """Update job status"""
    try:
        update_data = {'status': status, 'updated_at': 'now()'}
        if total_molecules is not None:
            update_data['total_molecules'] = total_molecules
        if completed_at:
            update_data['completed_at'] = completed_at
        supabase.table('jobs').update(update_data).eq('id', job_id).execute()
    except Exception as e:
        print(f"Error updating job {job_id}: {e}")
        raise

def get_job_by_id(job_id: str, user_id: str):
    """Get job by ID for specific user"""
    try:
        response = supabase.table('jobs').select('*').eq('id', job_id).eq('user_id', user_id).execute()
        if response.data:
            job = response.data[0]
            return {
                'id': job['id'],
                'user_id': job['user_id'],
                'parameters': job['parameters'],
                'status': job['status'],
                'total_molecules': job['total_molecules'],
                'created_at': job['created_at'],
                'updated_at': job['updated_at'],
                'completed_at': job['completed_at']
            }
        return None
    except Exception as e:
        print(f"Error getting job {job_id}: {e}")
        raise

def fetch_user_jobs(user_id: str, limit: int = 20, offset: int = 0):
    """Get user's jobs with pagination"""
    try:
        # Get total count
        count_response = supabase.table('jobs').select('id', count='exact').eq('user_id', user_id).execute()
        total_count = count_response.count

        # Get paginated jobs
        response = supabase.table('jobs').select('*').eq('user_id', user_id).order('created_at', desc=True).range(offset, offset + limit - 1).execute()

        jobs = response.data or []

        return {
            'jobs': [{
                'id': job['id'],
                'user_id': job['user_id'],
                'parameters': job['parameters'],
                'status': job['status'],
                'total_molecules': job['total_molecules'],
                'created_at': job['created_at'],
                'updated_at': job['updated_at'],
                'completed_at': job['completed_at']
            } for job in jobs],
            'total_count': total_count,
            'count': len(jobs),
            'limit': limit,
            'offset': offset,
            'has_more': (offset + limit) < total_count
        }
    except Exception as e:
        print(f"Error fetching user jobs: {e}")
        raise

def get_user_activity(user_id: str, limit: int = 10):
    """Get user's recent activity"""
    try:
        response = supabase.table('user_activity_history').select('*').eq('user_id', user_id).order('activity_date', desc=True).limit(limit).execute()
        activities = response.data or []
        return [{
            'activity_type': activity['activity_type'],
            'details': activity['details'],
            'credits_amount': activity['credits_amount'],
            'activity_date': activity['activity_date']
        } for activity in activities]
    except Exception as e:
        print(f"Error getting user activity: {e}")
        raise

def get_credit_history(user_id: str, limit: int = 10):
    """Get user's credit history"""
    try:
        response = supabase.table('credit_history').select('*').eq('user_id', user_id).order('created_at', desc=True).limit(limit).execute()
        history = response.data or []
        return [{
            'id': item['id'],
            'amount': item['amount'],
            'reason': item['reason'],
            'description': item['description'],
            'created_at': item['created_at']
        } for item in history]
    except Exception as e:
        print(f"Error getting credit history: {e}")
        raise

def add_credit_history(user_id: str, amount: int, reason: str, description: str):
    """Add credit history entry"""
    try:
        supabase.table('credit_history').insert({
            'user_id': user_id,
            'amount': amount,
            'reason': reason,
            'description': description
        }).execute()
    except Exception as e:
        print(f"Error adding credit history: {e}")
        raise


async def run_generation_with_timeout(params: dict) -> List[str]:
    loop = asyncio.get_running_loop()
    worker_fn = functools.partial(_run_generation_task, params)
    future = loop.run_in_executor(PROCESS_POOL, worker_fn)
    return await asyncio.wait_for(future, timeout=GENERATION_TIMEOUT_SECONDS)
def add_download_record(user_id: str, job_id: str, molecules_downloaded: int, credits_used: int, download_type: str):
    """Add download record"""
    try:
        supabase.table('downloads').insert({
            'user_id': user_id,
            'job_id': job_id,
            'molecules_downloaded': molecules_downloaded,
            'credits_used': credits_used,
            'download_type': download_type
        }).execute()
    except Exception as e:
        print(f"Error adding download record: {e}")
        raise

def add_user_activity(user_id: str, activity_type: str, details: str = None, credits_amount: int = None):
    """Add user activity entry (deprecated - using views now)"""
    # Since we now use views for activity history, this function is simplified
    # The activity is automatically tracked through jobs, downloads, and credit_history tables
    pass


def _run_generation_task(params: Dict[str, Any]) -> List[str]:
    """Isolated worker-safe wrapper for molecule generation."""
    return generate_functionalized_isomers(
        n_carbons=params.get("carbon_count", 0),
        functional_groups=params.get("functional_groups", []),
        n_double_bonds=params.get("double_bonds", 0),
        n_triple_bonds=params.get("triple_bonds", 0),
        n_rings=params.get("rings", 0),
        carbon_types=params.get("carbon_types") or ["primary", "secondary", "tertiary"],
    )


def cleanup_generation_memory():
    """Best-effort release of Python heap after heavy workloads."""
    gc.collect()

def store_job_results(user_id: str, job_id: str, molecules: List[str]):
    if not molecules:
        return
    try:
        payload = {"molecules": molecules}
        supabase.table('job_results').upsert({
            'job_id': job_id,
            'user_id': user_id,
            'payload': payload
        }).execute()
    except Exception as exc:
        print(f"Error storing job results for {job_id}: {exc}")


def load_job_results(job_id: str, user_id: str) -> List[str]:
    try:
        response = supabase.table('job_results').select('payload').eq('job_id', job_id).eq('user_id', user_id).execute()
        records = response.data or []
        if not records:
            return []
        payload = records[0].get('payload')
        if isinstance(payload, str):
            payload = json.loads(payload)
        return payload.get('molecules', []) if payload else []
    except Exception as exc:
        print(f"Error loading job results {job_id}: {exc}")
        return []


def delete_job_results(job_id: str):
    try:
        supabase.table('job_results').delete().eq('job_id', job_id).execute()
    except Exception as exc:
        print(f"Error deleting job results {job_id}: {exc}")


def prune_user_jobs(user_id: str, keep: int = 3):
    try:
        response = supabase.table('jobs').select('id').eq('user_id', user_id).order('created_at', desc=True).range(keep, keep + 99).execute()
        stale_jobs = response.data or []
        for job in stale_jobs:
            job_id = job['id']
            delete_job_results(job_id)
            supabase.table('jobs').delete().eq('id', job_id).execute()
    except Exception as exc:
        print(f"Error pruning jobs for user {user_id}: {exc}")

# Dependency to get current user from JWT token
def get_current_user(request: Request) -> str:
    """Extract user ID from Authorization header JWT token"""
    auth_header = request.headers.get('Authorization')
    if not auth_header or not auth_header.startswith('Bearer '):
        raise HTTPException(status_code=401, detail="Authorization header missing or invalid")

    token = auth_header.split(' ')[1]

    try:
        # Verify the JWT token with Supabase
        user_response = supabase.auth.get_user(token)
        return user_response.user.id
    except Exception as e:
        raise HTTPException(status_code=401, detail="Invalid token")

def get_user_tier(user_id: str) -> str:
    """Get user subscription tier"""
    user = get_user_by_id(user_id)
    if not user:
        return 'free'

    if user['is_fullaccess']:
        return 'fullaccess'
    elif user['credits'] >= 15:
        return 'paid'
    else:
        return 'free'


@app.post("/generate")
async def start_generation(request: JobRequest, user_id: str = Depends(get_current_user)):
    """Start molecule generation job (FREE)"""
    if STARTUP_ERROR:
        print("Startup error detected:\n", STARTUP_ERROR)
        raise HTTPException(status_code=500, detail="Server startup error: check application logs for import failures (RDKit/core).")
    # Check tier restrictions for large molecules
    user_tier = get_user_tier(user_id)
    if request.carbon_count > 7 and user_tier == 'free':
        raise HTTPException(
            status_code=403,
            detail="Generating molecules with more than 7 carbon atoms requires paid tier (15+ credits)"
        )

    # Server-side validation mirroring UI valency checks
    total_valency = 2 * request.carbon_count + 2 - 2 * request.double_bonds - 4 * request.triple_bonds
    if request.rings > 0:
        total_valency += request.rings

    functional_group_valency = len(request.functional_groups)  # Count of functional groups
    if functional_group_valency > total_valency:
        raise HTTPException(
            status_code=400,
            detail=f"Functional groups exceed available valency. Required: {functional_group_valency}, Available: {total_valency}"
        )

    # Validate structure possibility
    is_valid, error = validate_structure_possibility(
        request.carbon_count, request.functional_groups,
        request.double_bonds, request.triple_bonds, request.carbon_types, request.rings
    )
    if not is_valid:
        raise HTTPException(status_code=400, detail=error)

    # Create job in database
    job_id = str(uuid.uuid4())
    params_dict = request.model_dump()
    job_data = {
        'id': job_id,
        'user_id': user_id,
        'parameters': json.dumps(params_dict),
        'status': 'processing',
        'total_molecules': 0
    }

    try:
        create_job(job_data)
        add_user_activity(
            user_id,
            'job_started',
            f'Started generation: {request.carbon_count} carbons, {len(request.functional_groups)} functional groups',
            None
        )

        try:
            smiles_list = await run_generation_with_timeout(params_dict)
        except asyncio.TimeoutError:
            update_job_status(job_id, 'failed')
            raise HTTPException(status_code=504, detail="Generation timed out. Please simplify the request and try again.")
        except Exception as e:
            update_job_status(job_id, 'failed')
            raise HTTPException(status_code=500, detail=f"Generation failed: {str(e)}")

        completed_at = datetime.now().isoformat()
        update_job_status(job_id, 'completed', len(smiles_list), completed_at)

        store_job_results(user_id, job_id, smiles_list)

        prune_user_jobs(user_id)
        add_user_activity(user_id, 'job_completed', f'Generated {len(smiles_list)} molecules', None)

        return {
            "job_id": job_id,
            "status": "completed",
            "total_molecules": len(smiles_list),
            "molecules": smiles_list,
            "message": "Generation completed"
        }
    finally:
        cleanup_generation_memory()

@app.get("/jobs/{job_id}")
async def get_job_status(job_id: str, user_id: str = Depends(get_current_user)):
    """Get job status and basic info"""
    try:
        job = get_job_by_id(job_id, user_id)

        if not job:
            raise HTTPException(status_code=404, detail="Job not found")

        params = json.loads(job['parameters']) if job['parameters'] else {}
        molecules = load_job_results(job_id, user_id) if job['status'] == 'completed' else None

        return {
            "job_id": job['id'],
            "status": job['status'],
            "total_molecules": job['total_molecules'],
            "created_at": job['created_at'],
            "completed_at": job['completed_at'],
            "parameters": params,
            "molecules": molecules
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to get job: {str(e)}")

@app.get("/jobs/{job_id}/preview")
async def get_job_preview(job_id: str, user_id: str = Depends(get_current_user)):
    """Get preview of generated molecules (regenerate first few since not stored in DB)"""
    if STARTUP_ERROR:
        print("Startup error detected:\n", STARTUP_ERROR)
        raise HTTPException(status_code=500, detail="Server startup error: check application logs for import failures (RDKit/core).")
    try:
        # Get job
        job = get_job_by_id(job_id, user_id)

        if not job:
            raise HTTPException(status_code=404, detail="Job not found")

        if job['status'] != 'completed':
            raise HTTPException(status_code=400, detail="Job not completed yet")

        params = json.loads(job['parameters']) if job['parameters'] else {}
        smiles_list = load_job_results(job_id, user_id)

        if not smiles_list:
            smiles_list = generate_functionalized_isomers(
                n_carbons=params.get("carbon_count", 0),
                functional_groups=params.get("functional_groups", []),
                n_double_bonds=params.get("double_bonds", 0),
                n_triple_bonds=params.get("triple_bonds", 0),
                n_rings=params.get("rings", 0),
                carbon_types=params.get("carbon_types") or ["primary", "secondary", "tertiary"]
            )

        # Return first 3 molecules as preview
        preview_molecules = smiles_list[:3] if len(smiles_list) >= 3 else smiles_list

        return {
            "job_id": job_id,
            "status": "completed",
            "preview_molecules": preview_molecules,
            "preview_count": len(preview_molecules),
            "total_molecules": job['total_molecules'],
            "message": "Preview available - use download endpoint to get full results"
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to get preview: {str(e)}")

@app.post("/download")
async def download_molecules(request: DownloadRequest, user_id: str = Depends(get_current_user)):
    """Download molecules using credits (1 credit = 1000 molecules)"""
    if STARTUP_ERROR:
        print("Startup error detected:\n", STARTUP_ERROR)
        raise HTTPException(status_code=500, detail="Server startup error: check application logs for import failures (RDKit/core).")
    try:
        # Get job details
        job = get_job_by_id(request.job_id, user_id)

        if not job:
            raise HTTPException(status_code=404, detail="Job not found")

        if job['status'] != 'completed':
            raise HTTPException(status_code=400, detail="Job not completed")

        if request.molecules_count > job['total_molecules']:
            raise HTTPException(status_code=400, detail=f"Requested {request.molecules_count} molecules but only {job['total_molecules']} available")

        # Check download format permissions
        user_tier = get_user_tier(user_id)
        if request.download_format == 'all' and user_tier != 'fullaccess':
            raise HTTPException(status_code=403, detail="Full batch download requires fullaccess tier")

        # Calculate credits required (1 credit = 1000 molecules)
        credits_required = (request.molecules_count + 999) // 1000  # Ceiling division

        # Get user credits
        user_data = get_user_by_id(user_id)
        user_credits = user_data['credits']

        if user_credits < credits_required:
            raise HTTPException(
                status_code=402,
                detail=f"Insufficient credits. Required: {credits_required}, Available: {user_credits}"
            )

        # Deduct credits and add activity
        new_credits = user_credits - credits_required
        update_user_credits(user_id, new_credits)

        # Add download record
        download_format_for_record = 'sdf' if request.download_format == 'molsdf' else request.download_format
        add_download_record(user_id, request.job_id, request.molecules_count, credits_required, download_format_for_record)

        # Add credit history and activity
        add_credit_history(user_id, -credits_required, 'download', f"Downloaded {request.molecules_count} molecules in {request.download_format} format")
        add_user_activity(user_id, 'download',
                         f'Downloaded {request.molecules_count} molecules in {request.download_format} format',
                         -credits_required)

        params = json.loads(job['parameters']) if job['parameters'] else {}
        smiles_list = load_job_results(request.job_id, user_id)

        if not smiles_list:
            smiles_list = generate_functionalized_isomers(
                n_carbons=params.get("carbon_count", 0),
                functional_groups=params.get("functional_groups", []),
                n_double_bonds=params.get("double_bonds", 0),
                n_triple_bonds=params.get("triple_bonds", 0),
                n_rings=params.get("rings", 0),
                carbon_types=params.get("carbon_types") or ["primary", "secondary", "tertiary"]
            )

        # Take requested count
        selected_smiles = smiles_list[:request.molecules_count]

        # Format data based on download format
        if request.download_format == 'csv':
            data = 'SMILES\n' + '\n'.join(selected_smiles)
            content_type = 'text/csv'
            filename = f'molecules_{request.job_id[:8]}.csv'
        elif request.download_format == 'molsdf':
            # Create ZIP file with MOL, SDF, and CSV files
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                # Create CSV file
                csv_data = 'SMILES\n' + '\n'.join(selected_smiles)
                zip_file.writestr(f'molecules_{request.job_id[:8]}.csv', csv_data)
                
                # Create SDF file
                sdf_data = ''
                for i, smiles in enumerate(selected_smiles):
                    try:
                        mol = Chem.MolFromSmiles(smiles)
                        if mol:
                            mol_block = Chem.MolToMolBlock(mol)
                            sdf_data += mol_block + '\n$$$$\n'
                    except:
                        continue  # Skip invalid molecules
                zip_file.writestr(f'molecules_{request.job_id[:8]}.sdf', sdf_data)
                
                # Create individual MOL files
                for i, smiles in enumerate(selected_smiles):
                    try:
                        mol = Chem.MolFromSmiles(smiles)
                        if mol:
                            mol_block = Chem.MolToMolBlock(mol)
                            zip_file.writestr(f'molecule_{i+1:04d}.mol', mol_block)
                    except:
                        continue  # Skip invalid molecules
            
            data = base64.b64encode(zip_buffer.getvalue()).decode('utf-8')
            content_type = 'application/zip'
            filename = f'molecules_{request.job_id[:8]}_package.zip'
        else:  # json or all
            data = json.dumps({'molecules': selected_smiles})
            content_type = 'application/json'
            filename = f'molecules_{request.job_id[:8]}.json'

        return {
            "success": True,
            "credits_used": credits_required,
            "remaining_credits": new_credits,
            "data": data,
            "content_type": content_type,
            "filename": filename,
            "message": f"Downloaded {len(selected_smiles)} molecules. {credits_required} credits deducted."
        }

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Download failed: {str(e)}")

@app.get("/profile")
async def get_user_profile(user_id: str = Depends(get_current_user)):
    """Get user profile with credits and tier info"""
    try:
        # Get user data
        user_data = get_user_by_id(user_id)
        if not user_data:
            raise HTTPException(status_code=404, detail="User not found")

        # Get tier
        tier = get_user_tier(user_id)

        # Get recent activity (last 10 items)
        recent_activity = get_user_activity(user_id, 10)

        # Get credit history (last 10 items)
        credit_history = get_credit_history(user_id, 10)

        return {
            "user_id": user_data['id'],
            "email": user_data['email'],
            "credits": user_data['credits'],
            "subscription_tier": tier,
            "is_fullaccess": user_data['is_fullaccess'],
            "created_at": user_data['created_at'],
            "recent_activity": recent_activity,
            "credit_history": credit_history
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to get profile: {str(e)}")

@app.get("/jobs")
async def get_user_jobs(limit: int = 20, offset: int = 0, user_id: str = Depends(get_current_user)):
    """Get user's jobs history with pagination"""
    try:
        result = fetch_user_jobs(user_id, limit, offset)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to get jobs: {str(e)}")

@app.get("/jobs/summary")
async def get_jobs_summary(user_id: str = Depends(get_current_user)):
    """Get jobs summary for polling (serverless-friendly)"""
    try:
        # Get status counts using Supabase RPC or aggregation
        # Since Supabase doesn't have direct aggregation in Python client, we'll fetch and count
        response = supabase.table('jobs').select('status').eq('user_id', user_id).execute()
        jobs = response.data or []

        status_counts = {}
        for job in jobs:
            status = job['status']
            status_counts[status] = status_counts.get(status, 0) + 1

        # Ensure all status types are present
        for status in ['pending', 'processing', 'completed', 'failed']:
            if status not in status_counts:
                status_counts[status] = 0

        # Get active jobs (pending/processing)
        active_response = supabase.table('jobs').select('id,status,created_at').eq('user_id', user_id).in_('status', ['pending', 'processing']).order('created_at', desc=True).limit(5).execute()
        active_jobs = active_response.data or []

        active_jobs_formatted = [{
            'id': job['id'],
            'status': job['status'],
            'created_at': job['created_at']
        } for job in active_jobs]

        # Get total jobs count
        total_recent_jobs = len(jobs)

        return {
            "status_counts": status_counts,
            "active_jobs": active_jobs_formatted,
            "total_recent_jobs": total_recent_jobs
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to get jobs summary: {str(e)}")

@app.post("/credits/refill")
async def refill_credits(request: CreditRefillRequest, user_id: str = Depends(get_current_user)):
    """Admin endpoint to refill user credits (in production, integrate with payment system)"""
    try:
        # Get current credits
        user_data = get_user_by_id(user_id)
        current_credits = user_data['credits']

        # Add credits
        new_credits = current_credits + request.amount
        update_user_credits(user_id, new_credits)

        # Add credit history and activity
        add_credit_history(user_id, request.amount, 'refill', request.description)
        add_user_activity(user_id, 'credit_refill',
                         f'Added {request.amount} credits: {request.description}', request.amount)

        return {
            "success": True,
            "credits_added": request.amount,
            "new_balance": new_credits,
            "message": f"Added {request.amount} credits to account"
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to refill credits: {str(e)}")

@app.get("/", response_class=HTMLResponse)
async def root_status():
    return """
    <html>
        <head><title>Infinity Backend</title></head>
        <body>
            <h1>Infinity API is running</h1>
            <p>Status: OK</p>
        </body>
    </html>
    """
if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
