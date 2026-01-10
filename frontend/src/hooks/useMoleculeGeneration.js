import { useState, useEffect } from 'react'
import { API_BASE } from '../constants.js'
import { validateStructurePossibility } from '../utils/moleculeUtils.js'

export function useMoleculeGeneration(user, supabase, showToast) {
  const [formData, setFormData] = useState({
    carbon_count: 0,
    double_bonds: 0,
    triple_bonds: 0,
    rings: 0,
    carbon_types: ['primary', 'secondary', 'tertiary'],
    functional_groups: {}
  })
  const [currentJob, setCurrentJob] = useState(null)
  const [jobs, setJobs] = useState([])
  const [generating, setGenerating] = useState(false)

  const normalizeJob = (job = {}) => {
    const parsedParams = (() => {
      if (!job.parameters) return null
      if (typeof job.parameters === 'string') {
        try {
          return JSON.parse(job.parameters)
        } catch (_e) {
          return null
        }
      }
      return job.parameters
    })()

    return {
      job_id: job.job_id || job.id || null,
      status: job.status || 'pending',
      total_molecules: job.total_molecules || 0,
      created_at: job.created_at || null,
      completed_at: job.completed_at || null,
      parameters: parsedParams
    }
  }

  const handleInputChange = (field, value) => {
    setFormData(prev => ({ ...prev, [field]: value }))
  }

  const handleFunctionalGroupChange = (fg, value) => {
    const count = parseInt(value) || 0
    const carbonCount = formData.carbon_count
    const nDoubleBonds = formData.double_bonds
    const nTripleBonds = formData.triple_bonds
    const nRings = formData.rings

    let maxValency = 2 * carbonCount + 2 - 2 * nDoubleBonds - 4 * nTripleBonds
    if (nRings > 0) {
      maxValency += nRings
    }

    let currentValency = 0
    Object.entries(formData.functional_groups).forEach(([key, val]) => {
      if (key !== fg) currentValency += val
    })
    currentValency += count

    if (currentValency > maxValency) {
      showToast(`Cannot add more functional groups. The maximum allowed is ${maxValency}.`, 'error')
      return
    }

    setFormData(prev => ({
      ...prev,
      functional_groups: { ...prev.functional_groups, [fg]: count }
    }))
  }

  const handleCarbonTypeChange = (type, checked) => {
    setFormData(prev => ({
      ...prev,
      carbon_types: checked
        ? [...prev.carbon_types, type]
        : prev.carbon_types.filter(t => t !== type)
    }))
  }

  const generateMolecules = async () => {
    if (!user) {
      showToast('Please log in to generate molecules')
      return
    }

    const flattenedFunctionalGroups = Object.entries(formData.functional_groups)
      .filter(([_, count]) => count > 0)
      .flatMap(([fg, count]) => Array(count).fill(fg))

    const validation = validateStructurePossibility({
      carbonCount: formData.carbon_count,
      functionalGroups: flattenedFunctionalGroups,
      doubleBonds: formData.double_bonds,
      tripleBonds: formData.triple_bonds,
      rings: formData.rings,
      carbonTypes: formData.carbon_types
    })

    if (!validation.valid) {
      showToast(validation.message || 'Invalid molecular configuration', 'error')
      return
    }

    setGenerating(true)
    const payload = {
      carbon_count: formData.carbon_count,
      double_bonds: formData.double_bonds,
      triple_bonds: formData.triple_bonds,
      rings: formData.rings,
      carbon_types: formData.carbon_types,
      functional_groups: flattenedFunctionalGroups
    }

    try {
      const token = (await supabase.auth.getSession()).data.session?.access_token
      const response = await fetch(`${API_BASE}/generate`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${token}`
        },
        body: JSON.stringify(payload)
      })

      if (!response.ok) {
        const error = await response.json()
        throw new Error(error.detail || 'Generation failed')
      }

      const data = await response.json()
      setCurrentJob(normalizeJob({ job_id: data.job_id, status: 'pending', total_molecules: 0 }))
      showToast('Generation started! Check status below.', 'success')
      pollJobStatus(data.job_id)
      loadUserJobs()
    } catch (error) {
      showToast(error.message || 'Generation failed', 'error')
    } finally {
      setGenerating(false)
    }
  }

  const pollJobStatus = async (jobId) => {
    const interval = setInterval(async () => {
      try {
        const token = (await supabase.auth.getSession()).data.session?.access_token
        const response = await fetch(`${API_BASE}/jobs/${jobId}`, {
          headers: { 'Authorization': `Bearer ${token}` }
        })

        if (!response.ok) {
          clearInterval(interval)
          return
        }

        const data = await response.json()
        setCurrentJob(normalizeJob(data))

        if (data.status === 'completed' || data.status === 'failed') {
          clearInterval(interval)
          loadUserJobs()
          if (data.status === 'completed') {
            showToast(`Generation completed! ${data.total_molecules} molecules generated.`, 'success')
          } else if (data.status === 'failed') {
            showToast('Generation failed. Please review your inputs and try again.', 'error')
          }
        }
      } catch (error) {
        console.error('Polling error:', error)
        clearInterval(interval)
      }
    }, 2000)
  }

  const loadUserJobs = async () => {
    if (!user) return

    try {
      const token = (await supabase.auth.getSession()).data.session?.access_token
      const response = await fetch(`${API_BASE}/jobs`, {
        headers: { 'Authorization': `Bearer ${token}` }
      })

      if (response.ok) {
        const data = await response.json()
        const normalizedJobs = (data.jobs || []).map(normalizeJob)
        setJobs(normalizedJobs)

        if (!currentJob || ['completed', 'failed'].includes(currentJob.status)) {
          const latestCompleted = normalizedJobs.find(job => job.status === 'completed')
          if (latestCompleted) {
            setCurrentJob(latestCompleted)
          }
        }
      }
    } catch (error) {
      console.error('Failed to load jobs:', error)
    }
  }

  useEffect(() => {
    if (user) {
      loadUserJobs()
    }
  }, [user])

  return {
    formData,
    currentJob,
    jobs,
    generating,
    handleInputChange,
    handleFunctionalGroupChange,
    handleCarbonTypeChange,
    generateMolecules,
    loadUserJobs
  }
}