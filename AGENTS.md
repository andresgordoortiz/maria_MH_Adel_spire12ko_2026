# Environment isolation

- Use uv, pixi or containers for third-party tools and dependencies. Use a
  container for mixed Python/R/C++ software such as RNAmotifs.
- Do not install evaluation tools or dependencies into host system/user
  libraries or the existing betAS/analysis environments without explicit user
  approval. A temporary directory alone is not sufficient isolation.
- Review third-party setup/startup code for automatic dependency installation
  before execution. Keep source/reference mounts read-only and give containers
  write access only to their designated outputs.
