import { useState, useEffect } from 'react'
import { createClient } from '@supabase/supabase-js'
import './App.css'
import AuthModal from './components/AuthModal.jsx'
import ProfileModal from './components/ProfileModal.jsx'
import MoleculeForm from './components/MoleculeForm.jsx'

import DownloadSection from './components/DownloadSection.jsx'
import { useMoleculeGeneration } from './hooks/useMoleculeGeneration.js'
import { useUserProfile } from './hooks/useUserProfile.js'
import { FUNCTIONAL_GROUPS, API_BASE } from './constants.js'
import { updateFunctionalGroupStates, updateValencyStatus } from './utils/moleculeUtils.js'

let supabase
try {
  supabase = createClient(import.meta.env.VITE_SUPABASE_URL || 'dummy', import.meta.env.VITE_SUPABASE_ANON_KEY || 'dummy')
} catch (error) {
  console.error('Supabase init failed:', error)
  supabase = null
}

function App() {
  const [user, setUser] = useState(null)
  const [loading, setLoading] = useState(true)
  const [showAuthModal, setShowAuthModal] = useState(false)
  const [showProfile, setShowProfile] = useState(false)
  const [disabledGroups, setDisabledGroups] = useState({})
  const [valencyStatus, setValencyStatus] = useState('Valency status: OK')
  const [toast, setToast] = useState(null)
  const [downloading, setDownloading] = useState(false)

  const { profile, refillCredits } = useUserProfile(user, supabase)
  const {
    formData,
    currentJob,
    jobs,
    generating,
    handleInputChange,
    handleFunctionalGroupChange,
    handleCarbonTypeChange,
    generateMolecules
  } = useMoleculeGeneration(user, supabase, (message, type) => {
    setToast({ message, type })
    setTimeout(() => setToast(null), 5000)
  })

  const jobStatusMessage = (() => {
    if (!currentJob) return ''
    const shortId = currentJob.job_id ? currentJob.job_id.slice(0, 8) : ''
    switch (currentJob.status) {
      case 'pending':
        return `Job ${shortId} is queued. Preparing generation...`
      case 'processing':
        return `Job ${shortId} is generating molecules. Large batches can take a minute to render.`
      case 'completed':
        return `Job ${shortId} completed with ${currentJob.total_molecules} molecules. Ready to download.`
      case 'failed':
        return `Job ${shortId} failed. Please adjust inputs and try again.`
      default:
        return ''
    }
  })()

  useEffect(() => {
    const checkUser = async () => {
      if (supabase) {
        const { data: { user } } = await supabase.auth.getUser()
        setUser(user)
      }
      setLoading(false)
    }
    checkUser()

    let subscription = null
    if (supabase) {
      const { data: { subscription: sub } } = supabase.auth.onAuthStateChange(async (_event, session) => {
        setUser(session?.user ?? null)
        setLoading(false)
      })
      subscription = sub
    }

    return () => subscription?.unsubscribe()
  }, [])

  useEffect(() => {
    const disabled = updateFunctionalGroupStates(formData.carbon_count)
    setDisabledGroups(disabled)
    const status = updateValencyStatus(formData.carbon_count, formData.double_bonds, formData.triple_bonds, formData.rings, formData.functional_groups)
    setValencyStatus(status)
  }, [formData.carbon_count, formData.double_bonds, formData.triple_bonds, formData.rings, formData.functional_groups])

  const handleLogout = async () => {
    if (supabase) {
      await supabase.auth.signOut()
    }
  }

  const handleDownload = async (downloadForm) => {
    if (!user) {
      setShowAuthModal(true)
      return
    }

    setDownloading(true)
    try {
      const token = (await supabase.auth.getSession()).data.session?.access_token
      const response = await fetch(`${API_BASE}/download`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${token}`
        },
        body: JSON.stringify({
          job_id: downloadForm.jobId,
          molecules_count: downloadForm.moleculesCount,
          download_format: downloadForm.format
        })
      })

      if (response.status === 402) {
        // Insufficient credits - redirect to profile
        setShowProfile(true)
        setToast({ message: 'Insufficient credits. Please add more credits to download.', type: 'error' })
        setTimeout(() => setToast(null), 5000)
        return
      }

      if (!response.ok) {
        const error = await response.json()
        throw new Error(error.detail || 'Download failed')
      }

      const data = await response.json()

      // Create blob and download
      let blobData = data.data
      if (data.content_type === 'application/zip') {
        // Decode base64 for zip files
        blobData = Uint8Array.from(atob(data.data), c => c.charCodeAt(0))
      }
      const blob = new Blob([blobData], { type: data.content_type })
      const url = URL.createObjectURL(blob)
      const a = document.createElement('a')
      a.href = url
      a.download = data.filename
      document.body.appendChild(a)
      a.click()
      document.body.removeChild(a)
      URL.revokeObjectURL(url)

      setToast({ message: data.message, type: 'success' })
      setTimeout(() => setToast(null), 5000)
    } catch (error) {
      setToast({ message: 'Download failed: ' + error.message, type: 'error' })
      setTimeout(() => setToast(null), 5000)
    } finally {
      setDownloading(false)
    }
  }

  if (loading) return <div className="app"><div className="container">Loading...</div></div>

  return (
    <div className="app">
      <div className="profile-button" onClick={() => user ? setShowProfile(!showProfile) : setShowAuthModal(true)}>
        {user ? '👤' : '🔐'}
      </div>

      <AuthModal
        isOpen={showAuthModal}
        onClose={() => setShowAuthModal(false)}
        supabase={supabase}
      />

      <ProfileModal
        isOpen={showProfile}
        onClose={() => setShowProfile(false)}
        user={user}
        profile={profile}
        jobs={jobs}
        supabase={supabase}
        onLogout={handleLogout}
        onRefillCredits={refillCredits}
      />

      <div className="container">
        <div className="header">
          <h1>CHEM-∞ : The Universal Molecular Generator</h1>
        </div>
        <div className="system-info">System: 4 CPU cores, 16GB RAM</div>

        <MoleculeForm
          formData={formData}
          onInputChange={handleInputChange}
          onCarbonTypeChange={handleCarbonTypeChange}
          functionalGroups={FUNCTIONAL_GROUPS}
          disabledGroups={disabledGroups}
          onFunctionalGroupChange={handleFunctionalGroupChange}
          valencyStatus={valencyStatus}
        />

        <div className="buttons">
          <button>Clear Cache</button>
          <button onClick={generateMolecules} disabled={generating || !user}>
            {generating ? 'Generating...' : 'GENERATE COMPOUND STRUCTURES'}
          </button>
          <button onClick={() => setShowProfile(true)}>Profile</button>
        </div>

        {jobStatusMessage && (
          <div className={`job-status job-status-${currentJob?.status || 'default'}`}>
            {jobStatusMessage}
          </div>
        )}



        <DownloadSection
          user={user}
          jobs={jobs}
          profile={profile}
          currentJob={currentJob}
          onDownload={handleDownload}
          onShowProfile={() => setShowProfile(true)}
        />

        <div className="info-text">
          Bury the Worries of Drawing the Molecules<br />
          A code by Dr. Nitin Sapre, Computational Chemist<br />
          Create Chemical Database with MOL files and SDF files using parallel processing.
        </div>
      </div>

      {toast && (
        <div className={`toast toast-${toast.type}`}>
          {toast.message}
        </div>
      )}
    </div>
  )
}

export default App
