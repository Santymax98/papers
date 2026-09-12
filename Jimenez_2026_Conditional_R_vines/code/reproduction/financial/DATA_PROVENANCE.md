# Financial data provenance

The financial illustration uses the public Kenneth R. French Data Library file:

```text
10_Industry_Portfolios_daily_CSV.zip
https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/10_Industry_Portfolios_daily_CSV.zip
```

The source archive is not redistributed in this reviewer package because its
redistribution terms were not independently verified. The file used for the
reported analysis had the following identity:

```text
SHA-256: b1a23edba335ddd37223c46d75a6aa07827170dd6fd04d9a27249602e1054e4f
Size: 933519 bytes
Recorded local file timestamp: 2026-09-05T00:20:35Z
```

The fixed analysis window is 2015-01-01 through 2024-12-31. It retains 5,032
complete daily observations for the ten value-weighted industry portfolios:
NoDur, Durbl, Manuf, Enrgy, HiTec, Telcm, Shops, Hlth, Utils, and Other.

The deterministic processing is implemented in `run_financial.jl`: parse the
official daily CSV, retain the stated dates and complete rows, discard values
at or below -0.99 or non-finite values, then transform each margin to
mid-rank pseudo-observations divided by `n + 1`. Verify the downloaded file's
SHA-256 before rerunning the application.
