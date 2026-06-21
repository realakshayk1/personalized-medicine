"""Import side-effect module: importing this registers every primitive.

Importing each primitive module triggers its @primitive decorator, which inserts
the primitive into PRIMITIVE_REGISTRY. Consumers (the orchestrator, the planner,
the sandbox) import this single module instead of listing each primitive, so new
primitives become available everywhere by adding one line here.

Imports are alphabetical (import-sorter friendly). Each exists purely for its
registration side effect.
"""

import lattice_primitives.annotate.celltypist  # noqa: F401
import lattice_primitives.annotate.mllm_consensus  # noqa: F401
import lattice_primitives.cluster.leiden  # noqa: F401
import lattice_primitives.de.pseudobulk_deseq2  # noqa: F401
import lattice_primitives.de.rank_genes_groups  # noqa: F401
import lattice_primitives.dim_reduce.highly_variable_genes  # noqa: F401
import lattice_primitives.dim_reduce.neighbors  # noqa: F401
import lattice_primitives.dim_reduce.pca  # noqa: F401
import lattice_primitives.dim_reduce.umap  # noqa: F401
import lattice_primitives.integrate.harmony  # noqa: F401
import lattice_primitives.plot.dotplot  # noqa: F401
import lattice_primitives.plot.heatmap  # noqa: F401
import lattice_primitives.plot.plot_umap  # noqa: F401
import lattice_primitives.plot.volcano  # noqa: F401
import lattice_primitives.preprocess.filter_cells_basic  # noqa: F401
import lattice_primitives.preprocess.normalize_total_log1p  # noqa: F401
import lattice_primitives.preprocess.pseudobulk  # noqa: F401
import lattice_primitives.qc.calculate_qc_metrics  # noqa: F401
