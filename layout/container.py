import dash_bootstrap_components as dbc
from dash import dcc

from .plots import plots
from .sidebar import sidebar
from .analysis import main_analysis_tabs, side_analysis_tabs


main_body = dbc.Row(
    [
        dbc.Col(
            sidebar,
            className="p-0 sidebar pe-2",
            sm=12, md=12, lg=3
        ),
        dbc.Col(
            [
                plots,
                dbc.Row(
                    [
                        dbc.Col(
                            main_analysis_tabs,
                            # width=6,
                            className="block-display",
                            id="main-analysis-col"
                        ),
                        dbc.Col(
                            side_analysis_tabs,
                            # width=6,
                            className="no-display",
                            id="side-analysis-col"
                        )
                    ]
                )
            ],
            className="main-body mb-2",
            sm=12, md=12, lg=9
        ),
    ],
    className="g-0",
)

main_page = dbc.Col(
    [
        # Here we store any session variables that we need

        # Which plot is active. Value: 1 or 2 for main or side.
        dcc.Store(id="active-plot", storage_type="memory",
                  data=1, clear_data=True),
        # Whether we are in dual-mode or not. Value: 0 or 1.
        dcc.Store(id="dual-mode", storage_type="memory",
                  data=1, clear_data=True),

        # An output can only appear in a single dash callback. Therefore,
        # if there are multiple actions changing the same output, they all
        # need to be handled in a single function. To allow code readability
        # and extensibility, we split such actions into separate callbacks,
        # and set a "signal" output for each. These signals are then handled
        # by another callback that changes the desired output.

        # E.g., the dataset shape may change when 1) the data is loaded
        # 2) the plot is switched, 3) after preprocessing is run, etc.
        # For each such functionality we add a signal, and handle the
        # signals together in another callback.

        # Several UI components are duplicated (one for the main plot,
        # one for the side plot). Therefore, some signals are also duplicated,
        # one containing the prefix "main" and the other "side".

        # Signals data shape to change after loading data
        dcc.Store(id="shape-signal-load", storage_type="memory",
                  data=1, clear_data=True),
        # Signals data shape to change after running preprocessing
        dcc.Store(id="shape-signal-load-prep", storage_type="memory",
                  data=1, clear_data=True),
        # Signals data shape to change after running Bin To Gene conversion
        dcc.Store(id="shape-signal-load-prep-atac", storage_type="memory",
                  data=1, clear_data=True),
        # Signals data shape to change after switching plot
        dcc.Store(id="shape-signal-atoggle", storage_type="memory",
                  data=1, clear_data=True),

        # Used to clean analysis tables upon dataset switch
        dcc.Store(id="main-data-load-clean", storage_type="memory",
                  data=1, clear_data=True),
        dcc.Store(id="side-data-load-clean", storage_type="memory",
                  data=1, clear_data=True),

        # Signals cluster list to change after running any for of clustering,
        # including cluster, ss_cluster, label_transfer, subset_merge
        dcc.Store(id="main-cluster-list-signal", storage_type="memory",
                  data=1, clear_data=True),
        dcc.Store(id="side-cluster-list-signal", storage_type="memory",
                  data=1, clear_data=True),
        # Signals cluster list to change after merging a subset
        dcc.Store(id="main-cluster-merge-signal", storage_type="memory",
                  data=1, clear_data=True),
        dcc.Store(id="side-cluster-merge-signal", storage_type="memory",
                  data=1, clear_data=True),
        # Signals cluster list to change after adding a subset
        # (only applies to the DE cluster list)
        dcc.Store(id="main-subset-list-signal", storage_type="memory",
                  data=1, clear_data=True),
        dcc.Store(id="side-subset-list-signal", storage_type="memory",
                  data=1, clear_data=True),

        # Signal to update the plot, after merging a subset
        dcc.Store(id="merge-plot-signal", storage_type="memory",
                  data=1, clear_data=True),

        # Signal to update the feature list after loading data
        dcc.Store(id="feature-list-signal", storage_type="memory",
                  data=1, clear_data=True),
        # Signal to update the feature list after running preprocessing
        dcc.Store(id="feature-list-signal-prep", storage_type="memory",
                  data=1, clear_data=True),
        # Signal to update the feature list after running Bin To Gene
        dcc.Store(id="feature-list-signal-prep-atac", storage_type="memory",
                  data=1, clear_data=True),

        # Signal to add annotations to plot when hovering over cells
        dcc.Store(id="annotation-signal", storage_type="memory",
                  data=1, clear_data=True),

        # Signal to change main plot after running anything that changes it
        dcc.Store(id="main-plot-signal-code", storage_type="memory",
                  data=1, clear_data=True),
        dcc.Store(id="side-plot-signal-code", storage_type="memory",
                  data=1, clear_data=True),

        # Signal to change main plot after loading, running preprocessing,
        # or running Bin To Gene. If any embeddings / labels are found,
        # they will be displayed automatically.
        dcc.Store(id="data-loaded-plot-signal", storage_type="memory",
                  data=1, clear_data=True),
        dcc.Store(id="data-loaded-plot-signal-prep", storage_type="memory",
                  data=1, clear_data=True),
        dcc.Store(id="data-loaded-plot-signal-prep-atac",
                  storage_type="memory", data=1, clear_data=True),
        dcc.Store(id="data-loaded-annotation-table-signal",
                  storage_type="memory", data=1, clear_data=True),

        # Signal notification
        dcc.Store(id="push-notification",
                  storage_type="memory", data=1, clear_data=True),
        dcc.Store(id="main-other-features-change",
                  storage_type="memory", data=1, clear_data=True),
        dcc.Store(id="side-other-features-change",
                  storage_type="memory", data=1, clear_data=True),

        # Palette Signal
        dcc.Store(id="main-clear-palette-signal",
                  storage_type="memory", data=1, clear_data=True),
        dcc.Store(id="side-clear-palette-signal",
                  storage_type="memory", data=1, clear_data=True),
        dcc.Store(id="main-apply-palette-signal",
                  storage_type="memory", data=1, clear_data=True),
        dcc.Store(id="side-apply-palette-signal",
                  storage_type="memory", data=1, clear_data=True),

        dcc.Store(id="runtime-warning-signal",
                  storage_type="memory", data=1, clear_data=True),
        dcc.Store(id="shutdown-signal",
                  storage_type="memory", data=1, clear_data=True),
        dbc.Modal(
            [
                dbc.ModalHeader("Warning"),
                dbc.ModalBody(id="singler-runtime-body"),
                dbc.ModalFooter([
                    dbc.Button(
                        "Continue", id="singler-runtime-yes",
                        color="primary", className="ms-auto", n_clicks=0
                    ),
                    dbc.Button(
                        "Back", id="singler-runtime-no",
                        color="secondary", className="ms-auto", n_clicks=0
                    )
                ]),
            ],
            id="singler-runtime-modal",
            is_open=False
        ),
        # dcc.ConfirmDialog(id="singler-runtime-dialog"),
        main_body
    ],
    id = 'main-body',
    # className="no-display"
)

container=dbc.Row([main_page], id = 'container')
