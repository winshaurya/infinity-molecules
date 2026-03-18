const express = require('express');
const cors = require('cors');
const dotenv = require('dotenv');
const { createClient } = require('@supabase/supabase-js');
const { v4: uuidv4 } = require('uuid');

// Load environment variables
dotenv.config();

const app = express();

// CORS Configuration
const FRONTEND_URL = process.env.FRONTEND_URL || 'http://localhost:5173';
const BACKEND_URL = process.env.BACKEND_URL || 'http://localhost:8000';

app.use(cors({
  credentials: true,
  methods: ['GET', 'POST', 'PUT', 'DELETE', 'OPTIONS'],
  allowedHeaders: ['*']
}));

app.use(express.json());

// Generation runtime safeguards
const GENERATION_TIMEOUT_SECONDS = parseInt(process.env.GENERATION_TIMEOUT_SECONDS || '180');

// Supabase setup
const SUPABASE_URL = process.env.SUPABASE_URL;
const SUPABASE_SERVICE_ROLE_KEY = process.env.SUPABASE_SERVICE_ROLE_KEY;
const SUPABASE_ANON_KEY = process.env.SUPABASE_ANON_KEY;

console.log(`Loading SUPABASE_URL: ${SUPABASE_URL}`);
console.log(`Loading SUPABASE_SERVICE_ROLE_KEY: ${SUPABASE_SERVICE_ROLE_KEY ? SUPABASE_SERVICE_ROLE_KEY.substring(0, 20) + '...' : 'Not set'}`);
console.log(`Loading SUPABASE_ANON_KEY: ${SUPABASE_ANON_KEY ? SUPABASE_ANON_KEY.substring(0, 20) + '...' : 'Not set'}`);

const supabaseKey = SUPABASE_SERVICE_ROLE_KEY;

if (!SUPABASE_URL || !supabaseKey) {
  throw new Error('SUPABASE_URL and either SUPABASE_SERVICE_ROLE_KEY or SUPABASE_ANON_KEY must be set');
}

let supabase;
try {
  supabase = createClient(SUPABASE_URL, supabaseKey);
  console.log('✅ Supabase client initialized successfully');
} catch (e) {
  console.error(`❌ Supabase initialization failed: ${e}`);
  throw e;
}

// Import core logic (will implement later)
const {
  generateFunctionalizedIsomers,
  validateStructurePossibility,
  getFunctionalGroupType,
  getElementCounts,
  rdMolDescriptors
} = require('./core_logic');

// Database helper functions
async function getUserById(userId) {
  try {
    const { data, error } = await supabase
      .from('users')
      .select('*')
      .eq('id', userId)
      .single();

    if (error) throw error;

    if (data) {
      return {
        id: data.id,
        email: data.email,
        credits: data.credits,
        is_fullaccess: data.is_fullaccess,
        created_at: data.created_at,
        updated_at: data.updated_at
      };
    }
    return null;
  } catch (e) {
    console.error(`Error getting user ${userId}: ${e}`);
    throw e;
  }
}

async function updateUserCredits(userId, newCredits) {
  try {
    const { error } = await supabase
      .from('users')
      .update({ credits: newCredits, updated_at: 'now()' })
      .eq('id', userId);

    if (error) throw error;
  } catch (e) {
    console.error(`Error updating credits for user ${userId}: ${e}`);
    throw e;
  }
}

async function createJob(jobData) {
  try {
    const { error } = await supabase
      .from('jobs')
      .insert(jobData);

    if (error) throw error;
  } catch (e) {
    console.error(`Error creating job: ${e}`);
    throw e;
  }
}

async function updateJobStatus(jobId, status, totalMolecules = null, completedAt = null) {
  try {
    const updateData = { status, updated_at: 'now()' };
    if (totalMolecules !== null) updateData.total_molecules = totalMolecules;
    if (completedAt) updateData.completed_at = completedAt;

    const { error } = await supabase
      .from('jobs')
      .update(updateData)
      .eq('id', jobId);

    if (error) throw error;
  } catch (e) {
    console.error(`Error updating job ${jobId}: ${e}`);
    throw e;
  }
}

async function getJobById(jobId, userId) {
  try {
    const { data, error } = await supabase
      .from('jobs')
      .select('*')
      .eq('id', jobId)
      .eq('user_id', userId)
      .single();

    if (error) throw error;

    if (data) {
      return {
        id: data.id,
        user_id: data.user_id,
        parameters: data.parameters,
        status: data.status,
        total_molecules: data.total_molecules,
        created_at: data.created_at,
        updated_at: data.updated_at,
        completed_at: data.completed_at
      };
    }
    return null;
  } catch (e) {
    console.error(`Error getting job ${jobId}: ${e}`);
    throw e;
  }
}

async function fetchUserJobs(userId, limit = 20, offset = 0) {
  try {
    const { data: jobs, error: jobsError, count } = await supabase
      .from('jobs')
      .select('*', { count: 'exact' })
      .eq('user_id', userId)
      .order('created_at', { ascending: false })
      .range(offset, offset + limit - 1);

    if (jobsError) throw jobsError;

    return {
      jobs: jobs.map(job => ({
        id: job.id,
        user_id: job.user_id,
        parameters: job.parameters,
        status: job.status,
        total_molecules: job.total_molecules,
        created_at: job.created_at,
        updated_at: job.updated_at,
        completed_at: job.completed_at
      })),
      total_count: count,
      count: jobs.length,
      limit,
      offset,
      has_more: offset + limit < count
    };
  } catch (e) {
    console.error(`Error fetching user jobs: ${e}`);
    throw e;
  }
}

async function getUserActivity(userId, limit = 10) {
  try {
    const { data, error } = await supabase
      .from('user_activity_history')
      .select('*')
      .eq('user_id', userId)
      .order('activity_date', { ascending: false })
      .limit(limit);

    if (error) throw error;

    return data.map(activity => ({
      activity_type: activity.activity_type,
      details: activity.details,
      credits_amount: activity.credits_amount,
      activity_date: activity.activity_date
    }));
  } catch (e) {
    console.error(`Error getting user activity: ${e}`);
    throw e;
  }
}

async function getCreditHistory(userId, limit = 10) {
  try {
    const { data, error } = await supabase
      .from('credit_history')
      .select('*')
      .eq('user_id', userId)
      .order('created_at', { ascending: false })
      .limit(limit);

    if (error) throw error;

    return data.map(item => ({
      id: item.id,
      amount: item.amount,
      reason: item.reason,
      description: item.description,
      created_at: item.created_at
    }));
  } catch (e) {
    console.error(`Error getting credit history: ${e}`);
    throw e;
  }
}

async function addCreditHistory(userId, amount, reason, description) {
  try {
    const { error } = await supabase
      .from('credit_history')
      .insert({
        user_id: userId,
        amount,
        reason,
        description
      });

    if (error) throw error;
  } catch (e) {
    console.error(`Error adding credit history: ${e}`);
    throw e;
  }
}

async function addDownloadRecord(userId, jobId, moleculesDownloaded, creditsUsed, downloadType) {
  try {
    const { error } = await supabase
      .from('downloads')
      .insert({
        user_id: userId,
        job_id: jobId,
        molecules_downloaded: moleculesDownloaded,
        credits_used: creditsUsed,
        download_type: downloadType
      });

    if (error) throw error;
  } catch (e) {
    console.error(`Error adding download record: ${e}`);
    throw e;
  }
}

async function addUserActivity(userId, activityType, details = null, creditsAmount = null) {
  // Simplified - activity tracked through jobs, downloads, credit_history
}

async function storeJobResults(userId, jobId, molecules) {
  if (!molecules || molecules.length === 0) return;

  try {
    const payload = { molecules };
    const { error } = await supabase
      .from('job_results')
      .upsert({
        job_id: jobId,
        user_id: userId,
        payload
      });

    if (error) throw error;
  } catch (e) {
    console.error(`Error storing job results for ${jobId}: ${e}`);
  }
}

async function loadJobResults(jobId, userId) {
  try {
    const { data, error } = await supabase
      .from('job_results')
      .select('payload')
      .eq('job_id', jobId)
      .eq('user_id', userId)
      .single();

    if (error) throw error;

    if (data) {
      let payload = data.payload;
      if (typeof payload === 'string') {
        payload = JSON.parse(payload);
      }
      return payload.molecules || [];
    }
    return [];
  } catch (e) {
    console.error(`Error loading job results ${jobId}: ${e}`);
    return [];
  }
}

async function deleteJobResults(jobId) {
  try {
    const { error } = await supabase
      .from('job_results')
      .delete()
      .eq('job_id', jobId);

    if (error) throw error;
  } catch (e) {
    console.error(`Error deleting job results ${jobId}: ${e}`);
  }
}

async function pruneUserJobs(userId, keep = 3) {
  try {
    const { data: staleJobs, error } = await supabase
      .from('jobs')
      .select('id')
      .eq('user_id', userId)
      .order('created_at', { ascending: false })
      .range(keep, keep + 99);

    if (error) throw error;

    for (const job of staleJobs) {
      await deleteJobResults(job.id);
      await supabase.from('jobs').delete().eq('id', job.id);
    }
  } catch (e) {
    console.error(`Error pruning jobs for user ${userId}: ${e}`);
  }
}

// Middleware to get current user from JWT
async function getCurrentUser(req, res, next) {
  const authHeader = req.headers.authorization;
  if (!authHeader || !authHeader.startsWith('Bearer ')) {
    return res.status(401).json({ detail: 'Authorization header missing or invalid' });
  }

  const token = authHeader.split(' ')[1];

  try {
    const { data: { user }, error } = await supabase.auth.getUser(token);
    if (error) throw error;
    req.userId = user.id;
    next();
  } catch (e) {
    return res.status(401).json({ detail: 'Invalid token' });
  }
}

function getUserTier(userId) {
  // This will be implemented with caching later
  // For now, return 'free'
  return 'free';
}

const Joi = require('joi');

// Validation schemas
const jobRequestSchema = Joi.object({
  carbon_count: Joi.number().integer().min(0).max(20).required(),
  double_bonds: Joi.number().integer().min(0).max(10).required(),
  triple_bonds: Joi.number().integer().min(0).max(10).required(),
  rings: Joi.number().integer().min(0).max(10).required(),
  carbon_types: Joi.array().items(Joi.string().valid('primary', 'secondary', 'tertiary')).default(['primary', 'secondary', 'tertiary']),
  functional_groups: Joi.array().items(Joi.string()).max(50).required()
});

const downloadRequestSchema = Joi.object({
  job_id: Joi.string().required(),
  molecules_count: Joi.number().integer().min(1).required(),
  download_format: Joi.string().valid('csv', 'molsdf').required()
});

const creditRefillSchema = Joi.object({
  amount: Joi.number().integer().min(1).required(),
  description: Joi.string().default('Manual credit refill')
});

// Routes

// Generate endpoint
app.post('/generate', getCurrentUser, async (req, res) => {
  try {
    const { error, value } = jobRequestSchema.validate(req.body);
    if (error) {
      return res.status(400).json({ detail: error.details[0].message });
    }

    const params = value;
    const userId = req.userId;

    // Get user tier
    const userData = await getUserById(userId);
    let userTier = 'free';
    if (userData) {
      if (userData.is_fullaccess) {
        userTier = 'fullaccess';
      } else if (userData.credits >= 15) {
        userTier = 'paid';
      }
    }

    // Check tier restrictions
    if (params.carbon_count > 7 && userTier === 'free') {
      return res.status(403).json({
        detail: "Generating molecules with more than 7 carbon atoms requires paid tier (15+ credits)"
      });
    }

    // Server-side validation
    const totalValency = 2 * params.carbon_count + 2 - 2 * params.double_bonds - 4 * params.triple_bonds;
    const adjustedValency = params.rings > 0 ? totalValency + params.rings : totalValency;
    const functionalGroupValency = params.functional_groups.length;

    if (functionalGroupValency > adjustedValency) {
      return res.status(400).json({
        detail: `Functional groups exceed available valency. Required: ${functionalGroupValency}, Available: ${adjustedValency}`
      });
    }

    // Validate structure possibility
    const [isValid, errorMsg] = validateStructurePossibility(
      params.carbon_count, params.functional_groups,
      params.double_bonds, params.triple_bonds, params.carbon_types, params.rings
    );
    if (!isValid) {
      return res.status(400).json({ detail: errorMsg });
    }

    // Create job
    const jobId = uuidv4();
    const jobData = {
      id: jobId,
      user_id: userId,
      parameters: JSON.stringify(params),
      status: 'processing',
      total_molecules: 0
    };

    await createJob(jobData);

    // Generate molecules
    const smilesList = generateFunctionalizedIsomers(
      params.carbon_count,
      params.functional_groups,
      params.double_bonds,
      params.triple_bonds,
      params.rings,
      params.carbon_types
    );

    const completedAt = new Date().toISOString();
    await updateJobStatus(jobId, 'completed', smilesList.length, completedAt);

    await storeJobResults(userId, jobId, smilesList);
    await pruneUserJobs(userId);

    return res.json({
      job_id: jobId,
      status: 'completed',
      total_molecules: smilesList.length,
      molecules: smilesList,
      message: 'Generation completed'
    });

  } catch (e) {
    console.error('Generate error:', e);
    return res.status(500).json({ detail: `Generation failed: ${e.message}` });
  }
});

// Get job status
app.get('/jobs/:jobId', getCurrentUser, async (req, res) => {
  try {
    const { jobId } = req.params;
    const userId = req.userId;

    const job = await getJobById(jobId, userId);
    if (!job) {
      return res.status(404).json({ detail: 'Job not found' });
    }

    const params = JSON.parse(job.parameters || '{}');
    const molecules = job.status === 'completed' ? await loadJobResults(jobId, userId) : null;

    return res.json({
      job_id: job.id,
      status: job.status,
      total_molecules: job.total_molecules,
      created_at: job.created_at,
      completed_at: job.completed_at,
      parameters: params,
      molecules
    });

  } catch (e) {
    console.error('Get job error:', e);
    return res.status(500).json({ detail: `Failed to get job: ${e.message}` });
  }
});

// Get job preview
app.get('/jobs/:jobId/preview', getCurrentUser, async (req, res) => {
  try {
    const { jobId } = req.params;
    const userId = req.userId;

    const job = await getJobById(jobId, userId);
    if (!job) {
      return res.status(404).json({ detail: 'Job not found' });
    }

    if (job.status !== 'completed') {
      return res.status(400).json({ detail: 'Job not completed yet' });
    }

    const params = JSON.parse(job.parameters || '{}');
    const smilesList = await loadJobResults(jobId, userId);

    const previewMolecules = smilesList.slice(0, 3);

    return res.json({
      job_id: jobId,
      status: 'completed',
      preview_molecules: previewMolecules,
      preview_count: previewMolecules.length,
      total_molecules: job.total_molecules,
      message: 'Preview available'
    });

  } catch (e) {
    console.error('Get preview error:', e);
    return res.status(500).json({ detail: `Failed to get preview: ${e.message}` });
  }
});

// Download endpoint
app.post('/download', getCurrentUser, async (req, res) => {
  try {
    const { error, value } = downloadRequestSchema.validate(req.body);
    if (error) {
      return res.status(400).json({ detail: error.details[0].message });
    }

    const { job_id, molecules_count, download_format } = value;
    const userId = req.userId;

    const job = await getJobById(job_id, userId);
    if (!job) {
      return res.status(404).json({ detail: 'Job not found' });
    }

    if (job.status !== 'completed') {
      return res.status(400).json({ detail: 'Job not completed' });
    }

    if (molecules_count > job.total_molecules) {
      return res.status(400).json({
        detail: `Requested ${molecules_count} molecules but only ${job.total_molecules} available`
      });
    }

    const userData = await getUserById(userId);
    const creditsRequired = Math.ceil(molecules_count / 1000);

    if (userData.credits < creditsRequired) {
      return res.status(402).json({
        detail: `Insufficient credits. Required: ${creditsRequired}, Available: ${userData.credits}`
      });
    }

    // Deduct credits
    const newCredits = userData.credits - creditsRequired;
    await updateUserCredits(userId, newCredits);

    await addCreditHistory(userId, -creditsRequired, 'download', `Downloaded ${molecules_count} molecules in ${download_format} format`);
    await addDownloadRecord(userId, job_id, molecules_count, creditsRequired, download_format === 'molsdf' ? 'sdf' : download_format);

    const params = JSON.parse(job.parameters || '{}');
    const smilesList = await loadJobResults(job_id, userId);
    const selectedSmiles = smilesList.slice(0, molecules_count);

    if (download_format === 'csv') {
      const data = 'SMILES\n' + selectedSmiles.join('\n');
      return res.json({
        success: true,
        credits_used: creditsRequired,
        remaining_credits: newCredits,
        data: Buffer.from(data).toString('base64'),
        content_type: 'text/csv',
        filename: `molecules_${job_id.substring(0, 8)}.csv`,
        message: `Downloaded ${selectedSmiles.length} molecules. ${creditsRequired} credits deducted.`
      });
    } else if (download_format === 'molsdf') {
      // Simplified - return CSV for now
      const data = 'SMILES\n' + selectedSmiles.join('\n');
      return res.json({
        success: true,
        credits_used: creditsRequired,
        remaining_credits: newCredits,
        data: Buffer.from(data).toString('base64'),
        content_type: 'text/csv',
        filename: `molecules_${job_id.substring(0, 8)}.csv`,
        message: `Downloaded ${selectedSmiles.length} molecules. ${creditsRequired} credits deducted.`
      });
    }

  } catch (e) {
    console.error('Download error:', e);
    return res.status(500).json({ detail: `Download failed: ${e.message}` });
  }
});

// Get user profile
app.get('/profile', getCurrentUser, async (req, res) => {
  try {
    const userId = req.userId;
    const userData = await getUserById(userId);
    if (!userData) {
      return res.status(404).json({ detail: 'User not found' });
    }

    let tier = 'free';
    if (userData.is_fullaccess) {
      tier = 'fullaccess';
    } else if (userData.credits >= 15) {
      tier = 'paid';
    }

    const recentActivity = await getUserActivity(userId, 10);
    const creditHistory = await getCreditHistory(userId, 10);

    return res.json({
      user_id: userData.id,
      email: userData.email,
      credits: userData.credits,
      subscription_tier: tier,
      is_fullaccess: userData.is_fullaccess,
      created_at: userData.created_at,
      recent_activity: recentActivity,
      credit_history: creditHistory
    });

  } catch (e) {
    console.error('Profile error:', e);
    return res.status(500).json({ detail: `Failed to get profile: ${e.message}` });
  }
});

// Get user jobs
app.get('/jobs', getCurrentUser, async (req, res) => {
  try {
    const limit = parseInt(req.query.limit) || 20;
    const offset = parseInt(req.query.offset) || 0;
    const userId = req.userId;

    const result = await fetchUserJobs(userId, limit, offset);
    return res.json(result);

  } catch (e) {
    console.error('Get jobs error:', e);
    return res.status(500).json({ detail: `Failed to get jobs: ${e.message}` });
  }
});

// Get jobs summary
app.get('/jobs/summary', getCurrentUser, async (req, res) => {
  try {
    const userId = req.userId;

    // Simplified - return basic summary
    const result = await fetchUserJobs(userId, 100, 0);
    const statusCounts = { pending: 0, processing: 0, completed: 0, failed: 0 };

    for (const job of result.jobs) {
      statusCounts[job.status] = (statusCounts[job.status] || 0) + 1;
    }

    const activeJobs = result.jobs
      .filter(job => ['pending', 'processing'].includes(job.status))
      .slice(0, 5)
      .map(job => ({
        id: job.id,
        status: job.status,
        created_at: job.created_at
      }));

    return res.json({
      status_counts: statusCounts,
      active_jobs: activeJobs,
      total_recent_jobs: result.total_count
    });

  } catch (e) {
    console.error('Jobs summary error:', e);
    return res.status(500).json({ detail: `Failed to get jobs summary: ${e.message}` });
  }
});

// Refill credits
app.post('/credits/refill', getCurrentUser, async (req, res) => {
  try {
    const { error, value } = creditRefillSchema.validate(req.body);
    if (error) {
      return res.status(400).json({ detail: error.details[0].message });
    }

    const { amount, description } = value;
    const userId = req.userId;

    const userData = await getUserById(userId);
    const currentCredits = userData.credits;
    const newCredits = currentCredits + amount;

    await updateUserCredits(userId, newCredits);
    await addCreditHistory(userId, amount, 'refill', description);

    return res.json({
      success: true,
      credits_added: amount,
      new_balance: newCredits,
      message: `Added ${amount} credits to account`
    });

  } catch (e) {
    console.error('Refill credits error:', e);
    return res.status(500).json({ detail: `Failed to refill credits: ${e.message}` });
  }
});

// Root endpoint
app.get('/', (req, res) => {
  res.send(`
    <html>
        <head><title>Infinity Backend</title></head>
        <body>
            <h1>Infinity API is running</h1>
            <p>Status: OK</p>
        </body>
    </html>
  `);
});

const PORT = process.env.PORT || 8000;
app.listen(PORT, () => {
  console.log(`Server running on port ${PORT}`);
});

module.exports = app;