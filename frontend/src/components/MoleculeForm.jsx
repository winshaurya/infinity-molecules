function MoleculeForm({ formData, onInputChange, onCarbonTypeChange, functionalGroups, disabledGroups, onFunctionalGroupChange, valencyStatus }) {
  return (
    <>
      <div className="form-grid">
        <label>Number of Carbons (0-20):</label>
        <input
          type="number"
          min="0"
          max="20"
          value={formData.carbon_count}
          onChange={(e) => onInputChange('carbon_count', parseInt(e.target.value) || 0)}
        />

        <label>Number of Double Bonds:</label>
        <input
          type="number"
          min="0"
          max="10"
          value={formData.double_bonds}
          onChange={(e) => onInputChange('double_bonds', parseInt(e.target.value) || 0)}
        />

        <label>Number of Triple Bonds:</label>
        <input
          type="number"
          min="0"
          max="10"
          value={formData.triple_bonds}
          onChange={(e) => onInputChange('triple_bonds', parseInt(e.target.value) || 0)}
        />

        <label>Number of Rings:</label>
        <input
          type="number"
          min="0"
          max="10"
          value={formData.rings}
          onChange={(e) => onInputChange('rings', parseInt(e.target.value) || 0)}
        />

        <label>Carbon Types:</label>
        <div className="carbon-types">
          {['primary', 'secondary', 'tertiary'].map(type => (
            <label key={type}>
              <input
                type="checkbox"
                checked={formData.carbon_types.includes(type)}
                onChange={(e) => onCarbonTypeChange(type, e.target.checked)}
              />
              {type}
            </label>
          ))}
        </div>
      </div>

      <div className="functional-groups-title">Functional Groups (Counts):</div>
      <div className="functional-groups">
        {functionalGroups.map(fg => (
          <div key={fg} className="fg-item">
            <label>{fg}:</label>
            <input
              type="number"
              min="0"
              max="10"
              value={formData.functional_groups[fg] || 0}
              onChange={(e) => onFunctionalGroupChange(fg, e.target.value)}
              disabled={disabledGroups[fg]}
            />
          </div>
        ))}
      </div>

      <div className="valency-status">{valencyStatus}</div>

      <label>CPU Cores to use:</label>
      <div className="cpu-cores">
        <label><input type="radio" name="cpu" value="4" defaultChecked /> All 4 cores</label>
        <label><input type="radio" name="cpu" value="2" /> Half (2)</label>
        <label><input type="radio" name="cpu" value="custom" /> Custom:</label>
        <input type="number" min="1" max="4" defaultValue="4" style={{ width: '50px' }} />
      </div>
    </>
  )
}

export default MoleculeForm
