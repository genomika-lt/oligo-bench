import sys
import os
import shutil
import subprocess
import logging
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QPushButton,
    QProgressBar, QWidget, QMessageBox, QLabel, QHBoxLayout)
from PyQt6.QtCore import QThread, pyqtSignal, QTimer, Qt


class UninstallerThread(QThread):
    progress = pyqtSignal(int)
    log = pyqtSignal(str)
    error = pyqtSignal(str)

    def __init__(self, log_file,progress_bar):
        super().__init__()
        if getattr(sys, 'frozen', False):
            self.project_folder = os.path.dirname(sys.executable)
        else:
            self.project_folder = os.path.abspath(os.path.join(os.path.dirname(__file__)))
        self.project_env = "oligo"
        self.log_file = log_file
        self.progress_bar = progress_bar
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        handler = logging.FileHandler(log_file, mode='w')
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)

    def run(self):
        self._write_log("Checking if Conda environment 'test_env' exists.")
        self.progress_bar.setMaximum(100)
        result = subprocess.run(
            ["conda", "env", "list"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )

        if self.project_env in result.stdout:
            self._write_log("Conda environment 'test_env' found, proceeding with removal.")
            try:
                remove_result = subprocess.run(
                    ["conda", "env", "remove", "--prefix", os.path.join(self.project_folder, self.project_env), "-y"],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
                )

                if remove_result.returncode == 0:
                    self._write_log("Conda environment 'test_env' removed successfully.")
                else:
                    error_message = f"Failed to remove Conda environment: {remove_result.stderr}"
                    self._write_log(error_message)
                    self.error.emit(error_message)
            except subprocess.SubprocessError as e:
                error_message = f"Error removing Conda environment: {str(e)}"
                self._write_log(error_message)
                self.error.emit(error_message)
        else:
            self._write_log("Conda environment 'test_env' not found. Skipping removal.")

        self._write_log(f"Starting to remove files from: {self.project_folder}")
        total_items = sum([len(files) + len(dirs) for _, dirs, files in os.walk(self.project_folder)])
        removed_items = 0

        for root, dirs, files in os.walk(self.project_folder, topdown=False):
            for file in files:
                try:
                    file_path = os.path.join(root, file)
                    os.remove(file_path)
                    removed_items += 1
                    self.progress.emit(int((removed_items / total_items) * 50))
                except OSError as e:
                    self._write_log(f"Error removing file {file}: {str(e)}")
                    self.error.emit(f"Error removing file {file}: {str(e)}")

            for dir in dirs:
                try:
                    dir_path = os.path.join(root, dir)
                    shutil.rmtree(dir_path)
                    removed_items += 1
                    self.progress.emit(int((removed_items / total_items) * 50))
                except OSError as e:
                    self._write_log(f"Error removing directory {dir}: {str(e)}")
                    self.error.emit(f"Error removing directory {dir}: {str(e)}")

        try:
            if os.path.exists(self.project_folder):
                os.rmdir(self.project_folder)
                self._write_log("Project folder successfully deleted.")
            else:
                self._write_log("Project folder does not exist.")
        except OSError as e:
            self._write_log(f"Error deleting project folder: {str(e)}")
            self.error.emit(f"Error deleting project folder: {str(e)}")

        self.progress.emit(100)

    def _write_log(self, message):
        self.log.emit(message)
        self.logger.info(message)


class UninstallerApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Uninstaller App")
        self.setGeometry(300, 200, 600, 400)

        # Central widget and layout
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.layout = QVBoxLayout()
        self.central_widget.setLayout(self.layout)

        # Welcome label at the top
        self.welcome_label = QLabel("Welcome to Oligo-bench Uninstaller")
        self.welcome_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.layout.addWidget(self.welcome_label)

        # Progress bar setup
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)  # Initially invisible
        self.progress_bar.setValue(0)
        self.layout.addWidget(self.progress_bar)

        self.log_display = QLabel()
        self.layout.addWidget(self.log_display)

        button_layout = QHBoxLayout()

        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.clicked.connect(self.cancel_uninstall)
        self.cancel_button.hide()  # Initially hidden
        button_layout.addWidget(self.cancel_button)

        self.select_button = QPushButton("Start")
        self.select_button.clicked.connect(self.start_uninstallation)
        button_layout.addWidget(self.select_button)

        self.layout.addLayout(button_layout)

        self.uninstaller_thread = None
        self.log_file = "uninstaller.log"
        self.load_stylesheet()

    def start_uninstallation(self):
        self.uninstaller_thread = UninstallerThread(self.log_file,self.progress_bar)
        self.uninstaller_thread.progress.connect(self.update_progress)
        self.uninstaller_thread.log.connect(self.log_message)
        self.uninstaller_thread.error.connect(self.show_error)
        self.uninstaller_thread.finished.connect(self.on_uninstallation_finished)
        self.select_button.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.progress_bar.setMaximum(0)
        self.cancel_button.setVisible(True)
        self.uninstaller_thread.start()

    def cancel_uninstall(self):
        """
        Show a confirmation dialog to ask if the user really wants to cancel the uninstallation.
        If the user confirms, cancel the uninstallation and stop further progress.
        """
        reply = QMessageBox.question(
            self, "Cancel Uninstallation",
            "Are you sure you want to cancel the uninstallation?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No
        )

        if reply == QMessageBox.StandardButton.Yes:
            if self.uninstaller_thread:
                self.uninstaller_thread.terminate()
                self.log_message("Uninstallation process has been canceled.")
            self.cancel_button.hide()
            self.progress_bar.setVisible(False)
            self.select_button.setEnabled(True)
        else:
            return

    def load_stylesheet(self):
        """
        Loads the stylesheet from an external file.

        :return: None
        """
        if getattr(sys, 'frozen', False):
            css_path = os.path.join(sys._MEIPASS, "styles", "styles.qss")
        else:
            css_path = os.path.join(os.path.dirname(__file__), "styles",
                                    "styles.qss")
        try:
            with open(css_path, "r", encoding="utf-8") as file:
                stylesheet = file.read()
                self.setStyleSheet(stylesheet)
        except (FileNotFoundError, OSError) as e:
            self.show_error(f"Error loading stylesheet: {str(e)}")

    def update_progress(self, value):
        self.progress_bar.setValue(value)

    def log_message(self, message):
        self.log_display.setText(message)

    def show_error(self, error_message):
        QMessageBox.critical(self, "Error", error_message)
        self.log_display.setText(f"ERROR: {error_message}")

    def on_uninstallation_finished(self, timeout=3):
        """
        Show a notification message box with a dynamically changeable countdown.
        :param timeout: The number of seconds to wait before closing the message box.
        """
        msg_box = QMessageBox()
        msg_box.setIcon(QMessageBox.Icon.Information)
        msg_box.setWindowTitle("Uninstallation Complete")
        msg_box.setText(f"The project has been uninstalled successfully. Application will close in {timeout} seconds")
        msg_box.setStandardButtons(QMessageBox.StandardButton.Ok)

        def update_countdown():
            nonlocal timeout
            timeout -= 1
            msg_box.setText(f"The project has been uninstalled successfully. Application will close in {timeout} seconds")
            if timeout <= 0:
                msg_box.accept()
                countdown_timer.stop()
                QApplication.quit()

        countdown_timer = QTimer(self)
        countdown_timer.timeout.connect(update_countdown)
        countdown_timer.start(1000)
        msg_box.exec()

    def closeEvent(self, event):
        """
        Intercept the window close event to prevent the user from closing the app
        during the uninstallation process.
        """
        if self.uninstaller_thread and self.uninstaller_thread.isRunning():
            reply = QMessageBox.question(
                self, "Close Application",
                "Uninstallation is in progress. Are you sure you want to close the application?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No
            )

            if reply == QMessageBox.StandardButton.Yes:
                event.accept()
            else:
                event.ignore()
        else:
            event.accept()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = UninstallerApp()
    window.show()
    sys.exit(app.exec())
