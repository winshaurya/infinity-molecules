-- =====================================================
-- Chemistry SaaS Platform Database Schema
-- Supabase PostgreSQL
-- =====================================================
-- SUBSCRIPTION TIERS: free (<15 credits), paid (>=15 credits), fullaccess (special permission)
-- CREDITS: Used for DOWNLOADING molecules (1 credit = 1000 molecules), generation is FREE
-- HISTORY: Detailed tracking of generations, downloads, and credit refills
-- MOLECULES: NOT stored in database - only SMILES strings sent to frontend for rendering

-- Enable necessary extensions
CREATE EXTENSION IF NOT EXISTS "uuid-ossp";
CREATE EXTENSION IF NOT EXISTS "pgcrypto";

-- =====================================================
-- USERS TABLE
-- Stores user profile data and credits
-- =====================================================

CREATE TABLE IF NOT EXISTS public.users (
    id UUID REFERENCES auth.users(id) PRIMARY KEY,
    email TEXT NOT NULL,
    credits INTEGER DEFAULT 10 NOT NULL CHECK (credits >= 0),
    is_fullaccess BOOLEAN DEFAULT false NOT NULL,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT TIMEZONE('utc'::text, NOW()) NOT NULL,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT TIMEZONE('utc'::text, NOW()) NOT NULL
);
-- USERS TABLE: Stores user prcofiles with credits and fullaccess permission
-- Connects to: auth.users (Supabase's built-in authentication table)
-- Purpose: Manages user credits and special permissions
-- Subscription tiers calculated dynamically: free (<15 credits), paid (>=15 credits), fullaccess (flag)
-- Relationships: Referenced by jobs, downloads, credit_history tables
-- Constraints: Credits cannot be negative
-- is_fullaccess: Special permission for downloading entire batches

-- =====================================================
-- JOBS TABLE
-- Stores molecule generation jobs (FREE to generate)
-- =====================================================

CREATE TABLE IF NOT EXISTS public.jobs (
    id UUID DEFAULT uuid_generate_v4() PRIMARY KEY,
    user_id UUID REFERENCES public.users(id) ON DELETE CASCADE NOT NULL,
    parameters JSONB NOT NULL, -- Stores generation parameters
    status TEXT DEFAULT 'pending' CHECK (status IN ('pending', 'processing', 'completed', 'failed')),
    total_molecules INTEGER DEFAULT 0 CHECK (total_molecules >= 0),
    created_at TIMESTAMP WITH TIME ZONE DEFAULT TIMEZONE('utc'::text, NOW()) NOT NULL,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT TIMEZONE('utc'::text, NOW()) NOT NULL,
    completed_at TIMESTAMP WITH TIME ZONE
);
-- JOBS TABLE: Tracks molecule generation requests (FREE to generate)
-- Connects to: users table (user_id foreign key)
-- Purpose: Stores job parameters and generation status
-- Relationships: One user can have many jobs, referenced by downloads table
-- Status flow: pending -> processing -> completed (no more 'unlocked' status)
-- JSONB parameters: Stores carbons, bonds, functional groups, etc.
-- Generation is FREE: No credits required to generate molecules
-- Downloads cost credits: Tracked separately in downloads table

-- =====================================================
-- JOB RESULTS TABLE
-- Persists generated molecules per job (JSON payload)
-- =====================================================

CREATE TABLE IF NOT EXISTS public.job_results (
    job_id UUID PRIMARY KEY REFERENCES public.jobs(id) ON DELETE CASCADE,
    user_id UUID REFERENCES public.users(id) ON DELETE CASCADE NOT NULL,
    payload JSONB NOT NULL,
    result_size INTEGER GENERATED ALWAYS AS (jsonb_array_length(payload->'molecules')) STORED,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT TIMEZONE('utc'::text, NOW()) NOT NULL
);
-- JOB RESULTS TABLE: Stores generated SMILES payloads per job
-- Connects to: jobs table (job_id) and users table (user_id)
-- Purpose: Allows backend to fetch molecules without regenerating or relying on Storage
-- result_size: Derived column for quick counts and pruning heuristics
-- Lifecycle: Entries deleted when jobs are pruned to keep most recent generations only

-- =====================================================
-- DOWNLOADS TABLE
-- Tracks molecule downloads and credit usage
-- =====================================================

CREATE TABLE IF NOT EXISTS public.downloads (
    id UUID DEFAULT uuid_generate_v4() PRIMARY KEY,
    user_id UUID REFERENCES public.users(id) ON DELETE CASCADE NOT NULL,
    job_id UUID REFERENCES public.jobs(id) ON DELETE CASCADE NOT NULL,
    molecules_downloaded INTEGER NOT NULL CHECK (molecules_downloaded > 0),
    credits_used INTEGER NOT NULL CHECK (credits_used >= 0),
    download_type TEXT NOT NULL CHECK (download_type IN ('csv', 'sdf', 'all')),
    created_at TIMESTAMP WITH TIME ZONE DEFAULT TIMEZONE('utc'::text, NOW()) NOT NULL
);
-- DOWNLOADS TABLE: Detailed history of molecule downloads and credit usage
-- Connects to: users and jobs tables
-- Purpose: Tracks what users download from each generation job
-- molecules_downloaded: How many molecules were downloaded (up to total_molecules)
-- credits_used: Credits spent (1 credit = 1000 molecules)
-- download_type: 'csv', 'sdf', or 'all' (fullaccess users only)
-- History tracking: Shows user activity for profile/dashboard

-- =====================================================
-- CREDIT HISTORY TABLE
-- Tracks credit refills and usage
-- =====================================================

CREATE TABLE IF NOT EXISTS public.credit_history (
    id UUID DEFAULT uuid_generate_v4() PRIMARY KEY,
    user_id UUID REFERENCES public.users(id) ON DELETE CASCADE NOT NULL,
    amount INTEGER NOT NULL, -- Positive for refills, negative for usage
    reason TEXT NOT NULL CHECK (reason IN ('refill', 'download', 'admin')),
    description TEXT, -- Optional details
    created_at TIMESTAMP WITH TIME ZONE DEFAULT TIMEZONE('utc'::text, NOW()) NOT NULL
);
-- CREDIT HISTORY TABLE: Complete audit trail of credit changes
-- Connects to: users table
-- Purpose: Tracks all credit transactions for transparency
-- amount: Positive for credits added, negative for credits spent
-- reason: 'refill' (user added credits), 'download' (spent on downloads), 'admin' (manual adjustment)
-- description: Optional details about the transaction
-- Audit trail: Essential for user trust and debugging

-- =====================================================
-- INDEXES for Performance
-- =====================================================

-- Users table indexes
CREATE INDEX IF NOT EXISTS idx_users_email ON public.users(email);

-- Jobs table indexes
CREATE INDEX IF NOT EXISTS idx_jobs_user_id ON public.jobs(user_id);
CREATE INDEX IF NOT EXISTS idx_jobs_status ON public.jobs(status);
CREATE INDEX IF NOT EXISTS idx_jobs_created_at ON public.jobs(created_at DESC);

-- Job results indexes
CREATE INDEX IF NOT EXISTS idx_job_results_user_id ON public.job_results(user_id);
CREATE INDEX IF NOT EXISTS idx_job_results_created_at ON public.job_results(created_at DESC);

-- Downloads table indexes
CREATE INDEX IF NOT EXISTS idx_downloads_user_id ON public.downloads(user_id);
CREATE INDEX IF NOT EXISTS idx_downloads_job_id ON public.downloads(job_id);
CREATE INDEX IF NOT EXISTS idx_downloads_created_at ON public.downloads(created_at DESC);

-- Credit history indexes
CREATE INDEX IF NOT EXISTS idx_credit_history_user_id ON public.credit_history(user_id);
CREATE INDEX IF NOT EXISTS idx_credit_history_created_at ON public.credit_history(created_at DESC);
-- INDEXES: Optimize database query performance for common access patterns
-- User indexes: Fast email lookups and profile access
-- Job indexes: Speed up user job filtering, status queries, and chronological sorting
-- Download indexes: Track user download history and job-specific downloads
-- Credit history indexes: Audit trail queries and user credit transaction history
-- Performance impact: Critical for dashboard loading, history views, and real-time updates

-- =====================================================
-- ROW LEVEL SECURITY (RLS) POLICIES
-- Ensures users can only access their own data
-- =====================================================

-- Enable RLS on all tables
ALTER TABLE public.users ENABLE ROW LEVEL SECURITY;
ALTER TABLE public.jobs ENABLE ROW LEVEL SECURITY;
ALTER TABLE public.downloads ENABLE ROW LEVEL SECURITY;
ALTER TABLE public.credit_history ENABLE ROW LEVEL SECURITY;
ALTER TABLE public.job_results ENABLE ROW LEVEL SECURITY;

-- Users policies
CREATE POLICY "Users can view own profile" ON public.users
    FOR ALL USING (auth.uid() = id);

-- Jobs policies
CREATE POLICY "Users can view own jobs" ON public.jobs
    FOR ALL USING (auth.uid() = user_id);

-- Job results policies
CREATE POLICY "Users can manage own job results" ON public.job_results
    FOR ALL USING (auth.uid() = user_id);

-- Downloads policies
CREATE POLICY "Users can view own downloads" ON public.downloads
    FOR ALL USING (auth.uid() = user_id);

-- Credit history policies
CREATE POLICY "Users can view own credit history" ON public.credit_history
    FOR ALL USING (auth.uid() = user_id);
-- ROW LEVEL SECURITY: Ensures complete data isolation between users
-- ALTER TABLE...ENABLE RLS: Activates security policies on all tables
-- Users policy: Users can only access their own profile data
-- Jobs policy: Users can only see jobs they created
-- Downloads policy: Users can only see their download history
-- Credit history policy: Users can only see their credit transaction history
-- Security impact: Prevents users from accessing other users' data, even with direct database access
-- Supabase integration: auth.uid() returns the current authenticated user's ID from JWT token

-- =====================================================
-- FUNCTIONS and TRIGGERS
-- =====================================================

-- Function to update updated_at timestamp
CREATE OR REPLACE FUNCTION public.handle_updated_at()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = TIMEZONE('utc'::text, NOW());
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

-- Triggers for updated_at
CREATE TRIGGER handle_users_updated_at
    BEFORE UPDATE ON public.users
    FOR EACH ROW EXECUTE PROCEDURE public.handle_updated_at();

CREATE TRIGGER handle_jobs_updated_at
    BEFORE UPDATE ON public.jobs
    FOR EACH ROW EXECUTE PROCEDURE public.handle_updated_at();

-- Function to create user profile on signup
CREATE OR REPLACE FUNCTION public.handle_new_user()
RETURNS TRIGGER AS $$
BEGIN
    INSERT INTO public.users (id, email, credits)
    VALUES (NEW.id, NEW.email, 10);  -- Start with 10 credits
    RETURN NEW;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Trigger to create user profile on auth.users insert
CREATE TRIGGER on_auth_user_created
    AFTER INSERT ON auth.users
    FOR EACH ROW EXECUTE PROCEDURE public.handle_new_user();

-- Function to get user subscription tier
CREATE OR REPLACE FUNCTION public.get_user_tier(user_uuid UUID)
RETURNS TEXT AS $$
DECLARE
    user_credits INTEGER;
    is_fullaccess BOOLEAN;
BEGIN
    SELECT credits, is_fullaccess INTO user_credits, is_fullaccess
    FROM public.users WHERE id = user_uuid;

    IF is_fullaccess THEN
        RETURN 'fullaccess';
    ELSIF user_credits >= 15 THEN
        RETURN 'paid';
    ELSE
        RETURN 'free';
    END IF;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Function to download molecules and deduct credits
CREATE OR REPLACE FUNCTION public.download_molecules(
    user_uuid UUID,
    job_uuid UUID,
    download_count INTEGER,
    download_format TEXT
)
RETURNS JSON AS $$
DECLARE
    job_record RECORD;
    user_credits INTEGER;
    user_tier TEXT;
    credits_required INTEGER;
    download_id UUID;
    result JSON;
BEGIN
    -- Validate download format
    IF download_format NOT IN ('csv', 'sdf', 'all') THEN
        RETURN json_build_object('success', false, 'error', 'Invalid download format');
    END IF;

    -- Check if 'all' download requires fullaccess
    IF download_format = 'all' THEN
        user_tier := public.get_user_tier(user_uuid);
        IF user_tier != 'fullaccess' THEN
            RETURN json_build_object('success', false, 'error', 'Full batch download requires fullaccess tier');
        END IF;
    END IF;

    -- Get job details
    SELECT * INTO job_record FROM public.jobs WHERE id = job_uuid AND user_id = user_uuid;

    IF NOT FOUND THEN
        RETURN json_build_object('success', false, 'error', 'Job not found or access denied');
    END IF;

    IF job_record.status != 'completed' THEN
        RETURN json_build_object('success', false, 'error', 'Job must be completed to download');
    END IF;

    -- Validate download count doesn't exceed available molecules
    IF download_count > job_record.total_molecules THEN
        RETURN json_build_object('success', false, 'error',
                                format('Requested %s molecules but only %s available', download_count, job_record.total_molecules));
    END IF;

    -- Calculate credits required (1 credit = 1000 molecules)
    credits_required := CEIL(download_count::FLOAT / 1000);

    -- Get user credits
    SELECT credits INTO user_credits FROM public.users WHERE id = user_uuid;

    IF user_credits < credits_required THEN
        RETURN json_build_object('success', false, 'error', 'Insufficient credits',
                                'required', credits_required, 'available', user_credits);
    END IF;

    -- Generate download ID
    download_id := uuid_generate_v4();

    -- Atomic transaction: deduct credits and record download
    UPDATE public.users SET credits = credits - credits_required WHERE id = user_uuid;

    INSERT INTO public.downloads (id, user_id, job_id, molecules_downloaded, credits_used, download_type)
    VALUES (download_id, user_uuid, job_uuid, download_count, credits_required, download_format);

    INSERT INTO public.credit_history (user_id, amount, reason, description)
    VALUES (user_uuid, -credits_required, 'download',
            format('Downloaded %s molecules from job %s', download_count, job_uuid));

    -- Return success with download_id
    RETURN json_build_object('success', true, 'credits_used', credits_required,
                           'remaining_credits', user_credits - credits_required,
                           'download_id', download_id);

EXCEPTION
    WHEN OTHERS THEN
        RETURN json_build_object('success', false, 'error', SQLERRM);
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Function to refill user credits
CREATE OR REPLACE FUNCTION public.refill_credits(user_uuid UUID, credit_amount INTEGER, description TEXT DEFAULT NULL)
RETURNS JSON AS $$
DECLARE
    new_balance INTEGER;
BEGIN
    -- Update credits
    UPDATE public.users SET credits = credits + credit_amount WHERE id = user_uuid
    RETURNING credits INTO new_balance;

    -- Record in history
    INSERT INTO public.credit_history (user_id, amount, reason, description)
    VALUES (user_uuid, credit_amount, 'refill', description);

    RETURN json_build_object('success', true, 'credits_added', credit_amount, 'new_balance', new_balance);

EXCEPTION
    WHEN OTHERS THEN
        RETURN json_build_object('success', false, 'error', SQLERRM);
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;
-- FUNCTIONS & TRIGGERS: Automate common database operations
-- handle_updated_at(): PostgreSQL function that updates timestamp on row changes
-- BEFORE UPDATE triggers: Automatically set updated_at when users/jobs tables are modified
-- handle_new_user(): Creates user profile with 100 credits when new user signs up via Supabase Auth
-- AFTER INSERT trigger: Fires when new user is added to auth.users, creates corresponding profile
-- unlock_job_with_credits(): Atomic function to deduct credits and unlock job in single transaction
-- SECURITY DEFINER: Functions run with elevated privileges to modify user credits safely
-- Transaction safety: Credits deduction and job unlock happen atomically (all or nothing)
-- Error handling: Returns detailed JSON responses for different failure scenarios

-- =====================================================
-- VIEWS for Analytics and User Dashboard
-- =====================================================

-- View for user job and download statistics
CREATE OR REPLACE VIEW public.user_dashboard_stats AS
SELECT
    u.id as user_id,
    u.email,
    u.credits,
    public.get_user_tier(u.id) as subscription_tier,
    u.is_fullaccess,
    COUNT(DISTINCT j.id) as total_jobs,
    COUNT(DISTINCT CASE WHEN j.status = 'completed' THEN j.id END) as completed_jobs,
    SUM(j.total_molecules) as total_molecules_generated,
    COUNT(d.id) as total_downloads,
    SUM(d.molecules_downloaded) as total_molecules_downloaded,
    SUM(d.credits_used) as total_credits_spent_on_downloads,
    MAX(j.created_at) as last_job_date,
    MAX(d.created_at) as last_download_date
FROM public.users u
LEFT JOIN public.jobs j ON u.id = j.user_id
LEFT JOIN public.downloads d ON u.id = d.user_id
GROUP BY u.id, u.email, u.credits, u.is_fullaccess;

-- View for detailed user history (jobs + downloads + credit changes)
CREATE OR REPLACE VIEW public.user_activity_history AS
SELECT
    'job_created' as activity_type,
    j.user_id,
    j.id as reference_id,
    j.parameters::jsonb as details,
    j.created_at as activity_date,
    NULL as molecules_count,
    NULL as credits_amount
FROM public.jobs j

UNION ALL

SELECT
    'download' as activity_type,
    d.user_id,
    d.job_id as reference_id,
    json_build_object('molecules_downloaded', d.molecules_downloaded, 'download_type', d.download_type)::jsonb as details,
    d.created_at as activity_date,
    d.molecules_downloaded,
    -d.credits_used as credits_amount
FROM public.downloads d

UNION ALL

SELECT
    'credit_change' as activity_type,
    ch.user_id,
    ch.id as reference_id,
    json_build_object('reason', ch.reason, 'description', ch.description)::jsonb as details,
    ch.created_at as activity_date,
    NULL as molecules_count,
    ch.amount as credits_amount
FROM public.credit_history ch

ORDER BY activity_date DESC;
-- VIEWS: Pre-computed queries for analytics and user dashboard
-- user_dashboard_stats: Complete overview of user activity and current status
-- get_user_tier(): Dynamically calculates subscription tier (free/paid/fullaccess)
-- LEFT JOINs: Include users with no activity (shows 0 counts)
-- Aggregations: Sum molecules generated/downloaded, credits spent, etc.
-- user_activity_history: Unified timeline of all user activities
-- UNION ALL: Combines jobs, downloads, and credit changes into chronological order
-- Usage: Power user dashboards, analytics, and activity feeds

-- =====================================================
-- INITIAL DATA (Optional)
-- =====================================================

-- Insert default admin user (uncomment if needed)
-- INSERT INTO auth.users (id, email, encrypted_password, email_confirmed_at, created_at, updated_at)
-- VALUES (
--     'admin-uuid-here',
--     'admin@example.com',
--     crypt('admin_password', gen_salt('bf')),
--     NOW(),
--     NOW(),
--     NOW()
-- ) ON CONFLICT DO NOTHING;

-- =====================================================
-- GRANTS (if needed for specific roles)
-- =====================================================

-- Grant necessary permissions
GRANT USAGE ON SCHEMA public TO anon, authenticated;
GRANT ALL ON public.users TO anon, authenticated;
GRANT ALL ON public.jobs TO anon, authenticated;
GRANT ALL ON public.downloads TO anon, authenticated;
GRANT ALL ON public.credit_history TO anon, authenticated;
GRANT ALL ON public.user_dashboard_stats TO anon, authenticated;
GRANT ALL ON public.user_activity_history TO anon, authenticated;
-- GRANTS: Control access permissions for different user roles in Supabase
-- USAGE ON SCHEMA: Allows anon/authenticated users to access the public schema
-- GRANT ALL: Gives full CRUD permissions on tables to authenticated users
-- anon: Unauthenticated users (can read public data if RLS allows)
-- authenticated: Logged-in users (full access to their own data via RLS)
-- Includes all tables and views: users, jobs, downloads, credit_history, dashboard stats, activity history
-- Security: RLS policies still apply, so users only see their own data despite GRANT ALL

-- =====================================================
-- COMMENTS for Documentation
-- =====================================================

COMMENT ON TABLE public.users IS 'User profiles with credits and fullaccess permissions';
COMMENT ON TABLE public.jobs IS 'Molecule generation jobs (FREE to generate)';
COMMENT ON TABLE public.downloads IS 'Download history and credit usage tracking';
COMMENT ON TABLE public.credit_history IS 'Complete audit trail of credit transactions';
COMMENT ON VIEW public.user_dashboard_stats IS 'User dashboard statistics and subscription tier';
COMMENT ON VIEW public.user_activity_history IS 'Unified timeline of all user activities';

COMMENT ON COLUMN public.users.credits IS 'Available credits for downloading molecules (1 credit = 1000 molecules)';
COMMENT ON COLUMN public.users.is_fullaccess IS 'Special permission for full batch downloads';
COMMENT ON COLUMN public.jobs.parameters IS 'JSON object containing generation parameters (carbons, bonds, groups, etc.)';
COMMENT ON COLUMN public.jobs.status IS 'Job status: pending -> processing -> completed';
COMMENT ON COLUMN public.downloads.molecules_downloaded IS 'Number of molecules downloaded from this job';
COMMENT ON COLUMN public.downloads.credits_used IS 'Credits spent on this download (1 credit = 1000 molecules)';
COMMENT ON COLUMN public.downloads.download_type IS 'Download format: csv, sdf, or all (fullaccess only)';
COMMENT ON COLUMN public.credit_history.amount IS 'Positive for refills, negative for usage';
COMMENT ON COLUMN public.credit_history.reason IS 'refill, download, or admin';
-- COMMENTS: Human-readable documentation stored in PostgreSQL system catalogs
-- TABLE comments: Describe the purpose and contents of each table
-- VIEW comments: Explain what aggregated data the view provides
-- COLUMN comments: Clarify the meaning and format of specific columns
-- Usage: Comments appear in Supabase dashboard and can be queried via system tables
-- Maintenance: Comments help future developers understand the schema without reading code
-- Integration: These comments are used by ORMs and database tools for better UX
-- anon: Unauthenticated users (can read public data if RLS allows)
-- TABLE comments: Describe the purpose and contents of each table
-- anon: Unauthenticated users (can read public data if RLS allows)
