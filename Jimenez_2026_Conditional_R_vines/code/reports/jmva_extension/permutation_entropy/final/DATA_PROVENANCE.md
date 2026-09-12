# Data provenance

Source: Kenneth R. French Data Library, 10 Industry Portfolios, daily value-weighted returns.

- Source URL: https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/10_Industry_Portfolios_daily_CSV.zip
- Fixed local input: `data/10_Industry_Portfolios_daily_CSV.zip`
- Analysis dates: 2015-01-01 through 2024-12-31
- Retained observations: 5032
- First/last retained dates: 2015-01-02 / 2024-12-31
- Embedding: m=5, delay=1, overlapping windows
- Common marginal transform: mid-ranks over each complete series divided by n+1, followed by lagging
- Rounded ties and model fit: mid-ranks are used as a practical latent-continuous rank approximation, not as a claim of exact observed continuity
- Model selection: BIC is used as a working/composite pair-family score on the overlapping lag-vector rows; no iid likelihood-based selection inference is claimed
- Empirical ties: tied ordinal windows are excluded and reported in the summary; no iid inference is attached to overlapping frequencies
