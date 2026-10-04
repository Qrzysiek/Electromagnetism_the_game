# Audits

Each audit is archived here together with the commit it covered and the commit that
resolved it, so that the next audit can review only what changed since:
`git diff <fixed in>..HEAD`. Each folder holds the audit as delivered (paths inside its
documents refer to `audit/`, where it was written) and `RESOLUTION.md`, the verdict and
the fix for every remark.

| Audit | Audited commit | Fixed in | Folder | Next audit reviews |
|---|---|---|---|---|
| 2026-09-29, an external review (assisted by Grok 4.7) of the whole workspace | `25f821a` | `1292442` (2026-10-05) | [`2026-09-29_1292442/`](2026-09-29_1292442/RESOLUTION.md) | `git diff 1292442..HEAD` |
