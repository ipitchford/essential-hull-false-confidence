#!/bin/zsh
# Build paper.pdf from paper.md. Requires pandoc and tectonic. Run from the package root.
cd "${0:A:h}/.." && pandoc paper.md --from=markdown+tex_math_dollars+tex_math_single_backslash+autolink_bare_uris -s --pdf-engine=tectonic --include-in-header=header.tex -V geometry:margin=25mm -V fontsize=11pt -V colorlinks=true --resource-path=. -o paper.pdf
