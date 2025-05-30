import sys
from PySide6 import QtGui, QtWidgets
from PySide6.QtCore import QSize
from PySide6.QtGui import QIcon, QFont
from PySide6.QtWidgets import QApplication, QPushButton, QVBoxLayout, QWidget, QComboBox, QHBoxLayout, QListWidget, QRadioButton, QButtonGroup, QLabel, QCheckBox, QProgressBar
import pandas as pd
import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.figure import Figure


class Window(QWidget):

    def __init__(self, parent=None):
        super().__init__(parent)
        # First we set up all the buttons in the left-most columns
        gpp_folder_button = QPushButton(
            text="Select GPP data output directory", parent=self)
        gpp_folder_button.setFixedSize(200, 60)
        chip_file_button = QPushButton(text="Select corresponding .chip file",
                                       parent=self)
        chip_file_button.setFixedSize(200, 60)
        chip_file_colname_combobox = QComboBox()
        chip_file_colname_combobox.setFixedSize(200, 60)
        chip_file_colname_combobox.currentTextChanged.connect(
            self.chip_file_colname_combobox_text_changed)

        # Now we add some radio buttons for basic analysis options
        single_vs_matched_radio_group = QButtonGroup()
        single_radio = QRadioButton('single')
        single_vs_matched_radio_group.addButton(single_radio)
        matched_radio = QRadioButton('matched')
        single_vs_matched_radio_group.addButton(matched_radio)

        single_radio.clicked.connect(self.single_radio_clicked, )
        matched_radio.clicked.connect(self.matched_radio_clicked)

        radio_layout = QHBoxLayout()
        radio_layout.addWidget(QLabel("input type: "))
        radio_layout.addWidget(single_radio)
        radio_layout.addWidget(matched_radio)

        # Now a checkbox for Boolean option
        self.estimate_cells_checkbox = QCheckBox('estimate cells?')
        #estimate_cells.

        # Now a checkbox for Boolean option
        self.collapse_umis_checkbox = QCheckBox('collapse UMIs?')
        self.collapse_umis_checkbox.stateChanged.connect(self.load_gpp_data)
        # finally the run button
        run_button = QPushButton(text="Run FAUST", parent=self)
        run_button.setFixedSize(200, 60)
        plot_top_hits_button = QPushButton(text="Plot top hits", parent=self)
        plot_top_hits_button.setFixedSize(200, 60)

        input_list = QListWidget(parent=self)
        output_list = QListWidget(parent=self)
        comparisons_list = QListWidget(parent=self)
        self.input_list = input_list
        self.output_list = output_list
        self.comparisons_list = comparisons_list
        input_layout = QVBoxLayout()
        self.add_input_push_button = QPushButton(text="Add input", parent=self)
        self.add_input_push_button.clicked.connect(self.add_input)

        input_layout.addWidget(QLabel("Select input: "))
        input_layout.addWidget(self.input_list)
        input_layout.addWidget(self.add_input_push_button)
        output_layout = QVBoxLayout()

        self.add_output_push_button = QPushButton(text="Add output",
                                                  parent=self)
        self.add_output_push_button.clicked.connect(self.add_output)
        output_layout.addWidget(QLabel("Select output: "))

        output_layout.addWidget(self.output_list)
        output_layout.addWidget(self.add_output_push_button)

        comparisons_layout = QVBoxLayout()
        comparisons_layout.addWidget(QLabel("Comparisons: "))
        comparisons_layout.addWidget(self.comparisons_list)
        clear_button = QPushButton(text="Clear comparison list", parent=self)
        clear_button.clicked.connect(self.clear_comparisons)

        comparisons_layout.addWidget(clear_button)

        buttonlayout = QVBoxLayout()
        buttonlayout.addWidget(gpp_folder_button)
        buttonlayout.addWidget(chip_file_button)
        buttonlayout.addWidget(chip_file_colname_combobox)
        buttonlayout.addLayout(radio_layout)
        checkboxlayout = QHBoxLayout()

        checkboxlayout.addWidget(self.estimate_cells_checkbox)
        checkboxlayout.addWidget(self.collapse_umis_checkbox)
        buttonlayout.addLayout(checkboxlayout)
        buttonlayout.addWidget(run_button)
        buttonlayout.addWidget(plot_top_hits_button)

        # the spreadlayout is just the 4 columns
        spreadlayout = QHBoxLayout()
        spreadlayout.addLayout(buttonlayout)
        spreadlayout.addLayout(input_layout)
        spreadlayout.addLayout(output_layout)
        spreadlayout.addLayout(comparisons_layout)
        gpp_folder_button.clicked.connect(self.load_gpp_results)
        chip_file_button.clicked.connect(self.load_chip_file)
        run_button.clicked.connect(self.run_FAUST)
        plot_top_hits_button.clicked.connect(self.run_plot_top_hits)

        self.chip_file_colname_combobox = chip_file_colname_combobox
        self.gpp_folder_button = gpp_folder_button
        self.chip_file_button = chip_file_button
        self.run_button = run_button
        self.plot_top_hits_button = plot_top_hits_button
        self.folderpath = None
        self.chipfile = None
        self.chip_df = None
        self.chip_file_colname = None
        self.add_input_push_button.setDisabled(True)
        self.add_output_push_button.setDisabled(True)
        self.run_button.setDisabled(True)
        self.plot_top_hits_button.setDisabled(True)

        # comparison lists:
        self.input_type = None
        self.inputs = []
        self.outputs = []

        # targets that may be selected as controls
        self.targets = []
        self.controls = []
        self.target_list = QListWidget(parent=self)
        self.target_list.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        target_layout = QVBoxLayout()
        target_layout.addWidget(QLabel("Targets (select controls): "))
        target_layout.addWidget(self.target_list)
        select_targets_as_controls_button = QPushButton(
            text="Select targets as controls", parent=self)
        select_targets_as_controls_button.clicked.connect(self.add_controls)
        target_layout.addWidget(select_targets_as_controls_button)
        spreadlayout.addLayout(target_layout)
        # selected controls
        self.controls_list = QListWidget(parent=self)
        controls_layout = QVBoxLayout()
        controls_layout.addWidget(QLabel("Controls: "))
        controls_layout.addWidget(self.controls_list)
        clear_controls_button = QPushButton(text="Clear controls", parent=self)
        clear_controls_button.clicked.connect(self.clear_controls)
        controls_layout.addWidget(clear_controls_button)
        spreadlayout.addLayout(controls_layout)

        # samples
        self.replicate_groups = []
        self.samples_list = QListWidget(parent=self)
        self.samples_list.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        samples_to_group_layout = QVBoxLayout()
        samples_to_group_layout.addWidget(
            QLabel("Sample (group replicates together): "))
        samples_to_group_layout.addWidget(self.samples_list)
        group_samples_as_replicates_button = QPushButton(
            text="Group samples as replicates", parent=self)
        group_samples_as_replicates_button.clicked.connect(
            self.group_replicates)
        samples_to_group_layout.addWidget(group_samples_as_replicates_button)
        spreadlayout.addLayout(samples_to_group_layout)

        # replicate groups
        self.replicate_groups = []
        self.replicate_groups_list = QListWidget(parent=self)
        replicate_groups_layout = QVBoxLayout()
        replicate_groups_layout.addWidget(QLabel("Replicate groups: "))
        replicate_groups_layout.addWidget(self.replicate_groups_list)
        clear_replicate_groups_button = QPushButton(
            text="Clear replicate groups", parent=self)
        clear_replicate_groups_button.clicked.connect(
            self.clear_replicate_groups)
        replicate_groups_layout.addWidget(clear_replicate_groups_button)
        spreadlayout.addLayout(replicate_groups_layout)

        # add in progress bar
        fulllayout = QVBoxLayout()
        fulllayout.addLayout(spreadlayout)
        self.progress_bar = QProgressBar(self)
        fulllayout.addWidget(self.progress_bar)
        self.setLayout(fulllayout)

    def preflight_check(self):
        checks = [
            self.folderpath is None, self.chipfile is None,
            self.chip_file_colname is None,
            len(self.controls) == 0,
            len(self.inputs) == 0,
            len(self.outputs) == 0
        ]
        if np.any(checks):
            self.run_button.setDisabled(True)
        else:
            self.run_button.setDisabled(False)

    def load_gpp_data(self):

        if (self.folderpath is not None) and (self.chipfile is not None) and (
                self.chip_file_colname is not None):
            from faust.utilities import read_gpp_output
            collapse_umis = self.collapse_umis_checkbox.isChecked()
            gpp = read_gpp_output(
                self.folderpath,
                chipfile=self.chipfile,
                chipfile_gene_symbol_colname=self.chip_file_colname,
                collapse_umis=collapse_umis)
            self.gpp = gpp
            self.input_list.clear()
            self.output_list.clear()
            self.samples_list.clear()
            for col in gpp.columns:
                self.input_list.addItem(col)
                self.output_list.addItem(col)
                self.samples_list.addItem(col)
            self.preflight_check()

    def rename_and_color_button(self, button, newstring):
        max_length = int(0.1 * (button.size().width()) - 3)
        if len(newstring) > max_length:
            newstring_abbreviated = '...' + newstring[::-1][0:max_length][::-1]
        else:
            newstring_abbreviated = newstring
        button.setText(newstring_abbreviated)
        button.setStyleSheet("background-color : green")

    def load_gpp_results(self):
        dialog = QtWidgets.QFileDialog()
        folderpath = QtWidgets.QFileDialog.getExistingDirectory(
            self, 'Select Folder')
        self.folderpath = folderpath
        self.rename_and_color_button(self.gpp_folder_button, folderpath)

        self.load_gpp_data()
        self.preflight_check()

    def load_chip_file(self):
        dialog = QtWidgets.QFileDialog()
        fileName, _filter = QtWidgets.QFileDialog.getOpenFileName(
            self, "Select .chip file", ".", "(*.chip)")
        self.chipfile = fileName
        self.chip_df = pd.read_table(fileName)
        self.chip_file_colname_combobox.clear()
        for col in [
                x for x in self.chip_df.columns if x != 'Barcode Sequence'
        ]:
            self.chip_file_colname_combobox.addItem(col)

        self.rename_and_color_button(self.chip_file_button, fileName)

        self.load_gpp_data()
        self.preflight_check()

    def chip_file_colname_combobox_text_changed(self, s):
        self.chip_file_colname = s
        self.chip_file_colname_combobox.setStyleSheet(
            "background-color : green")
        self.populate_targets()
        self.load_gpp_data()
        self.preflight_check()

    def populate_targets(self):
        self.targets = list(np.unique(self.chip_df[self.chip_file_colname]))
        self.target_list.clear()
        for target in self.targets:
            self.target_list.addItem(target)

    def single_radio_clicked(self):
        self.add_input_push_button.setDisabled(False)
        self.add_output_push_button.setDisabled(False)
        self.add_output_push_button.setText('Add input')

        self.add_output_push_button.setText('Add output')
        self.input_type = 'single'

    def matched_radio_clicked(self):
        self.add_input_push_button.setDisabled(True)
        self.add_output_push_button.setDisabled(False)
        self.add_output_push_button.setText('Add matched input/output pair')
        self.input_type = 'matched'

    def add_input(self):
        for item in self.input_list.selectedItems():
            self.inputs.append(item.text())
        self.arrange_comparisons()

    def add_output(self):
        if self.input_type == 'single':
            for item in self.output_list.selectedItems():
                self.outputs.append(item.text())
            self.arrange_comparisons()
        elif self.input_type == 'matched':
            for item in self.input_list.selectedItems():
                self.inputs.append(item.text())
            for item in self.output_list.selectedItems():
                self.outputs.append(item.text())
            self.arrange_comparisons()

    def arrange_comparisons(self):
        if self.input_type == 'single':
            inputs = '+'.join(self.inputs)
            self.comparisons_list.clear()
            for output in self.outputs:
                comparison = '{} --> {}'.format(inputs, output)
                self.comparisons_list.addItem(comparison)
        elif self.input_type == 'matched':
            self.comparisons_list.clear()
            for i, o in zip(self.inputs, self.outputs):
                comparison = '{} --> {}'.format(i, o)
                self.comparisons_list.addItem(comparison)
        self.preflight_check()

    def clear_comparisons(self):
        self.inputs = []
        self.outputs = []
        self.arrange_comparisons()

    def add_controls(self):
        for item in self.target_list.selectedItems():
            self.controls.append(item.text())
        self.controls = list(np.unique(self.controls))
        for control in self.controls:
            self.controls_list.addItem(control)
        self.preflight_check()

    def clear_controls(self):
        self.controls_list.clear()

    def group_replicates(self):
        group = []
        for replicate in self.samples_list.selectedItems():
            group.append(replicate.text())
        self.replicate_groups.append(group)
        self.replicate_groups_list.clear()
        for replicate_group in self.replicate_groups:
            self.replicate_groups_list.addItem(', '.join(replicate_group))

    def clear_replicate_groups(self):
        self.replicate_groups = []
        self.replicate_groups_list.clear()

    def run_FAUST(self):
        from faust.utilities import get_summary_df
        output_filename, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save File",
            "",
            "Excel file or text table (*.xlsx *txt *csv *tsv);;All Files (*)",
        )
        if not np.any([
                output_filename.endswith(y)
                for y in ['.xlsx', '.txt', '.csv', '.tsv']
        ]):
            output_filename += '.txt'

        self.summary_df = get_summary_df(
            self.gpp,
            self.controls,
            self.inputs,
            self.outputs,
            input_type=self.input_type,
            estimate_cells=self.estimate_cells_checkbox.isChecked(),
            verbose=False,
            progress_logger=self.progress_bar)

        if len(self.replicate_groups) > 0:
            self.summary_df['input->output'] = self.summary_df.apply(
                lambda x: '{}->{}'.format(x['input_site'], x['output_site']),
                axis=1)
            comparison2equivalence_class = {}
            for comparison in self.summary_df['input->output'].unique():
                equivalence_class = comparison
                for group in self.replicate_groups:
                    for replicate in group:
                        if replicate in comparison:
                            equivalence_class = equivalence_class.replace(
                                replicate, '|'.join(group))
                comparison2equivalence_class[comparison] = equivalence_class
            self.summary_df['comparison_equivalence_class'] = self.summary_df[
                'input->output'].apply(lambda x: comparison2equivalence_class[
                    x] if x in comparison2equivalence_class.keys() else x)
            from faust.utilities import get_replicate_aggregated_statistics
            get_replicate_aggregated_statistics(
                self.summary_df,
                aggregation_column='comparison_equivalence_class',
                inplace=True)

        if output_filename.endswith('.xlsx'):
            self.summary_df.to_excel(output_filename, index=False)
        elif output_filename.endswith('.tsv'):
            self.summary_df.to_csv(output_filename, sep='\t', index=False)
        else:
            self.summary_df.to_csv(output_filename, index=False)
        self.plot_top_hits_button.setDisabled(False)

    def run_plot_top_hits(self):
        from faust.visualization import plot_top_hits
        import matplotlib
        matplotlib.use('QtAgg')
        #matplotlib.rcParams['backend.qt5']='PySide'
        fig, ax = plot_top_hits(self.summary_df, both_sides=True)
        fig.show()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec())
