# CHANGELOG

## [3.0.0] - 2025-07-09

### Major Changes
- Complete implementation of all 28 ACMG/AMP evidence criteria (PVS1, PS1-4, PM1-6, PP1-5, BA1, BS1-4, BP1-7)
- Replaced AI-based placeholder logic with interactive user evaluation system
- Refactored codebase with modular architecture under `src/` directory
- Added comprehensive test suite for all criteria validation

### Added
- Interactive evaluation for criteria requiring literature review (PS1, PS4, PM3, PM5, PP4, BP2, BP5, BP6)
- Enhanced metascore calculation using multiple computational predictors
- Population frequency analysis with gene-specific thresholds
- SpliceAI integration for splice impact assessment
- Test mode for automated validation of all criteria
- Detailed evidence breakdown in classification reports

### Enhanced
- Input validation and error handling
- Report generation with clinical recommendations
- Documentation with complete criteria descriptions
- User interface for interactive criteria evaluation

### Fixed
- Report output path inconsistency between Python and executable versions
- Input validation for edge cases and missing data
- Memory management during large dataset processing
- Interactive workflow user experience

### Removed
- AI-based automatic literature scanning functionality
- Mock phenotype similarity scoring system
- Placeholder implementations for incomplete criteria

---

## [2.2.0] - Previous Release
### Added
- Enhanced de novo variant analysis for PS2/PM6 criteria
- Improved alias handling system for variant notation
- Population frequency analysis framework

### Fixed
- Input validation improvements
- Report generation error handling

---

## [2.1.0] - Previous Release
### Added
- Basic ACMG criteria implementation framework
- Computational metascore integration
- Interactive mode foundation

### Enhanced
- User interface improvements
- Data validation system
