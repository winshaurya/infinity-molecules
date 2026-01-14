import { useState, useEffect } from 'react'

function DownloadSection({ user, jobs, profile, onDownload, onRefillCredits, currentJob, onShowProfile }) {
  const [downloadForm, setDownloadForm] = useState({
    jobId: currentJob?.job_id || '',
    moleculesCount: currentJob?.total_molecules ? Math.min(1000, currentJob.total_molecules) : 1000,
    format: 'csv'
  })

  // Update jobId when currentJob changes
  useEffect(() => {
    setDownloadForm(prev => ({
      ...prev,
      jobId: currentJob?.job_id || '',
      moleculesCount: currentJob?.total_molecules
        ? Math.min(prev.moleculesCount || 1000, currentJob.total_molecules)
        : (prev.moleculesCount || 1000)
    }))
  }, [currentJob?.job_id, currentJob?.total_molecules])
  const [isDownloading, setIsDownloading] = useState(false)
  const [downloadProgress, setDownloadProgress] = useState(0)

  const maxMolecules = currentJob?.total_molecules || 0

  const handleDownload = async () => {
    setIsDownloading(true)
    setDownloadProgress(0)

    try {
      await onDownload(downloadForm, (progress) => {
        setDownloadProgress(progress)
      })
      setDownloadProgress(100)
    } finally {
      setTimeout(() => {
        setIsDownloading(false)
        setDownloadProgress(0)
      }, 1000)
    }
  }

  if (!user) return null

  return (
    <div className="download-section">
      <h3>Download Molecules</h3>
      <div className="download-form">
        {currentJob && currentJob.status === 'completed' ? (
          <div className="current-job-info">
            <p><strong>Current Job:</strong> {currentJob.job_id.slice(0, 8)} - {currentJob.total_molecules} molecules</p>
            <input type="hidden" value={currentJob.job_id} />
          </div>
        ) : (
          <p className="no-current-job">No completed job available for download</p>
        )}
        <input
          type="number"
          placeholder="Molecules to download"
          value={downloadForm.moleculesCount}
          onChange={(e) => {
            const rawValue = parseInt(e.target.value, 10) || 0
            setDownloadForm(prev => ({ ...prev, moleculesCount: Math.max(0, rawValue) }))
          }}
          min="0"
          step="1000"
        />
        <select
          value={downloadForm.format}
          onChange={(e) => setDownloadForm(prev => ({ ...prev, format: e.target.value }))}
        >
          <option value="molsdf">MOL/SDF Package (ZIP)</option>
          <option value="csv">CSV Only</option>
        </select>
        <button onClick={handleDownload} disabled={!(currentJob && currentJob.status === 'completed' && currentJob.total_molecules > 0) || isDownloading}>
          {isDownloading ? 'Downloading...' : 'Download'}
        </button>
        {downloadForm.moleculesCount > 0 && (
          <span className="credit-cost">
            Cost: {Math.ceil(downloadForm.moleculesCount / 1000)} credits
          </span>
        )}
      </div>

      {isDownloading && downloadForm.format === 'molsdf' && (
        <div className="download-progress">
          <div className="progress-text">
            Rendering compounds... {downloadProgress}%
          </div>
          <div className="progress-bar">
            <div
              className="progress-fill"
              style={{ width: `${downloadProgress}%` }}
            ></div>
          </div>
        </div>
      )}
    </div>
  )
}

export default DownloadSection
