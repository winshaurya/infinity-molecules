import { useState } from 'react'

function ProfileModal({ isOpen, onClose, user, profile, jobs, supabase, onLogout, onRefillCredits }) {
  const [profileTab, setProfileTab] = useState('profile')
  const [passwordData, setPasswordData] = useState({ current: '', new: '', confirm: '' })
  const [toast, setToast] = useState(null)

  const showToast = (message, type = 'error') => {
    setToast({ message, type })
    setTimeout(() => setToast(null), 5000)
  }

  const handleChangePassword = async (e) => {
    e.preventDefault()
    if (passwordData.new !== passwordData.confirm) {
      showToast('New passwords do not match')
      return
    }
    const { error } = await supabase.auth.updateUser({
      password: passwordData.new
    })
    if (error) {
      showToast('Password change failed: ' + error.message)
    } else {
      showToast('Password changed successfully!', 'success')
      setPasswordData({ current: '', new: '', confirm: '' })
    }
  }

  if (!isOpen) return null

  return (
    <div className="modal-overlay" onClick={onClose}>
      <div className="profile-modal" onClick={(e) => e.stopPropagation()}>
        <div className="profile-header">
          <h2>Account Settings</h2>
          <button onClick={onClose} className="close-btn">&times;</button>
        </div>
        <div className="profile-tabs">
          <button
            className={profileTab === 'profile' ? 'active' : ''}
            onClick={() => setProfileTab('profile')}
          >
            Profile
          </button>
          <button
            className={profileTab === 'history' ? 'active' : ''}
            onClick={() => setProfileTab('history')}
          >
            History
          </button>
          <button
            className={profileTab === 'password' ? 'active' : ''}
            onClick={() => setProfileTab('password')}
          >
            Password
          </button>
          <button
            className={profileTab === 'payment' ? 'active' : ''}
            onClick={() => setProfileTab('payment')}
          >
            Credits
          </button>
        </div>
        <div className="profile-content">
          {profileTab === 'profile' && (
            <div className="profile-section">
              <h3>Profile Information</h3>
              <div className="profile-field">
                <label>Email:</label>
                <span>{user?.email || 'Not logged in'}</span>
              </div>
              <div className="profile-field">
                <label>Credits:</label>
                <span>{profile?.credits || 0}</span>
              </div>
              <div className="profile-field">
                <label>Subscription Tier:</label>
                <span className={`tier-${profile?.subscription_tier || 'free'}`}>
                  {profile?.subscription_tier || 'free'}
                </span>
              </div>
              <div className="profile-field">
                <label>Full Access:</label>
                <span>{profile?.is_fullaccess ? 'Yes' : 'No'}</span>
              </div>
              <div className="profile-field">
                <label>Member Since:</label>
                <span>{profile?.created_at ? new Date(profile.created_at).toLocaleDateString() : 'N/A'}</span>
              </div>
              <button onClick={onLogout} className="logout-btn">Logout</button>
            </div>
          )}
          {profileTab === 'history' && (
            <div className="profile-section">
              <h3>Activity History</h3>
              <div className="activity-timeline">
                {profile?.recent_activity?.map((activity, idx) => (
                  <div key={idx} className="timeline-item">
                    <div className="timeline-marker">
                      <div className={`activity-icon ${activity.activity_type}`}>
                        {activity.activity_type === 'job_created' && '🧪'}
                        {activity.activity_type === 'job_completed' && '✅'}
                        {activity.activity_type === 'download' && '📥'}
                        {activity.activity_type === 'credit_change' && '💰'}
                        {!['job_created', 'job_completed', 'download', 'credit_change'].includes(activity.activity_type) && '📋'}
                      </div>
                    </div>
                    <div className="timeline-content">
                      <div className="activity-header">
                        <span className="activity-type-label">
                          {activity.activity_type.replace('_', ' ').replace(/\b\w/g, l => l.toUpperCase())}
                        </span>
                        <span className="activity-time">
                          {new Date(activity.activity_date).toLocaleDateString()} at {new Date(activity.activity_date).toLocaleTimeString([], {hour: '2-digit', minute:'2-digit'})}
                        </span>
                      </div>
                      <div className="activity-body">
                        {activity.activity_type === 'job_created' && activity.details ? (
                          <div className="job-params">
                            <div className="param-grid">
                              {(() => {
                                try {
                                  const params = typeof activity.details === 'string' ? JSON.parse(activity.details) : activity.details;
                                  return (
                                    <>
                                      <div className="param-item">
                                        <span className="param-label">Carbons:</span>
                                        <span className="param-value">{params.carbon_count}</span>
                                      </div>
                                      <div className="param-item">
                                        <span className="param-label">Double bonds:</span>
                                        <span className="param-value">{params.double_bonds}</span>
                                      </div>
                                      <div className="param-item">
                                        <span className="param-label">Triple bonds:</span>
                                        <span className="param-value">{params.triple_bonds}</span>
                                      </div>
                                      <div className="param-item">
                                        <span className="param-label">Rings:</span>
                                        <span className="param-value">{params.rings}</span>
                                      </div>
                                      <div className="param-item full-width">
                                        <span className="param-label">Functional groups:</span>
                                        <span className="param-value">{params.functional_groups?.join(', ') || 'None'}</span>
                                      </div>
                                    </>
                                  );
                                } catch (e) {
                                  return <span>{activity.details}</span>;
                                }
                              })()}
                            </div>
                          </div>
                        ) : activity.activity_type === 'download' && activity.details ? (
                          <div className="download-info">
                            {(() => {
                              try {
                                const details = typeof activity.details === 'string' ? JSON.parse(activity.details) : activity.details;
                                return (
                                  <span>
                                    Downloaded {details.molecules_downloaded} molecules in {details.download_type} format
                                  </span>
                                );
                              } catch (e) {
                                return <span>{activity.details}</span>;
                              }
                            })()}
                          </div>
                        ) : activity.activity_type === 'credit_change' && activity.details ? (
                          <div className="credit-info">
                            {(() => {
                              try {
                                const details = typeof activity.details === 'string' ? JSON.parse(activity.details) : activity.details;
                                return (
                                  <span>
                                    {details.reason === 'refill' ? 'Credits added' : 'Credits used'}: {details.description || 'No description'}
                                  </span>
                                );
                              } catch (e) {
                                return <span>{activity.details}</span>;
                              }
                            })()}
                          </div>
                        ) : (
                          <span>{activity.details}</span>
                        )}
                      </div>
                      {activity.credits_amount !== null && (
                        <div className={`activity-credits ${activity.credits_amount > 0 ? 'positive' : 'negative'}`}>
                          {activity.credits_amount > 0 ? '+' : ''}{activity.credits_amount} credits
                        </div>
                      )}
                    </div>
                  </div>
                )) || (
                  <div className="no-activity">
                    <div className="no-activity-icon">📭</div>
                    <p>No recent activity</p>
                    <small>Start generating molecules to see your activity here!</small>
                  </div>
                )}
              </div>
            </div>
          )}
          {profileTab === 'password' && (
            <div className="profile-section">
              <h3>Change Password</h3>
              <form onSubmit={handleChangePassword}>
                <div className="form-group">
                  <label>New Password:</label>
                  <input
                    type="password"
                    value={passwordData.new}
                    onChange={(e) => setPasswordData(prev => ({ ...prev, new: e.target.value }))}
                    required
                  />
                </div>
                <div className="form-group">
                  <label>Confirm New Password:</label>
                  <input
                    type="password"
                    value={passwordData.confirm}
                    onChange={(e) => setPasswordData(prev => ({ ...prev, confirm: e.target.value }))}
                    required
                  />
                </div>
                <button type="submit">Update Password</button>
              </form>
            </div>
          )}
          {profileTab === 'payment' && (
            <div className="profile-section">
              <h3>Credits & Billing</h3>
              <div className="profile-field">
                <label>Current Credits:</label>
                <span>{profile?.credits || 0}</span>
              </div>
              <div className="profile-field">
                <label>Subscription Tier:</label>
                <span className={`tier-${profile?.subscription_tier || 'free'}`}>
                  {profile?.subscription_tier || 'free'}
                </span>
              </div>
              <div className="credits-info">
                <p><strong>How credits work:</strong></p>
                <ul>
                  <li>1 credit = 1000 molecules downloaded</li>
                  <li>Free tier: less than 15 credits</li>
                  <li>Paid tier: &ge;15 credits (can generate till 7 carbon atoms only )</li>
                  <li>Full Access tier: Download entire batches</li>
                </ul>
              </div>
              <div className="credit-actions">
                <button onClick={async () => {
                  const result = await onRefillCredits(50)
                  if (result?.success) {
                    showToast(result.message, 'success')
                  } else {
                    showToast(result?.message || 'Failed to add credits', 'error')
                  }
                }} className="refill-btn">Add 50 Credits (Test)</button>
                <button onClick={async () => {
                  const result = await onRefillCredits(100)
                  if (result?.success) {
                    showToast(result.message, 'success')
                  } else {
                    showToast(result?.message || 'Failed to add credits', 'error')
                  }
                }} className="refill-btn">Add 100 Credits (Test)</button>
              </div>

              <div className="credit-history-section">
                <h4>Credit History</h4>
                <div className="credit-history-list">
                  {profile?.credit_history?.map((item, idx) => (
                    <div key={idx} className="credit-history-item">
                      <div className="credit-history-header">
                        <span className={`credit-amount ${item.amount > 0 ? 'positive' : 'negative'}`}>
                          {item.amount > 0 ? '+' : ''}{item.amount} credits
                        </span>
                        <span className="credit-reason">{item.reason}</span>
                        <span className="credit-date">
                          {new Date(item.created_at).toLocaleDateString()} at {new Date(item.created_at).toLocaleTimeString([], {hour: '2-digit', minute:'2-digit'})}
                        </span>
                      </div>
                      {item.description && (
                        <div className="credit-description">{item.description}</div>
                      )}
                    </div>
                  )) || <p>No credit history available</p>}
                </div>
              </div>
            </div>
          )}
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

export default ProfileModal