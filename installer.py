import os
import subprocess
import sys
import logging
import urllib
import zipfile
from pathlib import Path
import requests
from PyQt6.QtGui import QFont

from PyQt6.QtWidgets import QApplication, QWidget, QVBoxLayout, QPushButton, QProgressBar, QLabel, QFileDialog, \
    QHBoxLayout, QMessageBox, QMainWindow, QStackedWidget, QLineEdit
from PyQt6.QtCore import QThread, pyqtSignal, Qt

log_file="installation.log"
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.FileHandler(log_file)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
with open(log_file, 'w'):
    logger.debug(f"Log file {log_file} initialized.")

class DownloadThread(QThread):
    progress = pyqtSignal(int)
    status = pyqtSignal(str)
    complete = pyqtSignal(bool)
    indeterminate = pyqtSignal()

    def __init__(self, url, save_path, parent=None):
        super().__init__(parent)
        self.url = url
        self.save_path = Path(save_path)
        self.cancelled = False

    def run(self):
        """Main function to run the download, unzip, and installation process."""
        logger.info(f"Download started for URL: {self.url} to path: {self.save_path}")

        self.status.emit("Starting download...")
        self.save_path.mkdir(parents=True, exist_ok=True)
        zip_file_name = self.url.split("/")[-1]
        zip_file_path = self.save_path / zip_file_name
        self.download_file(self.url,zip_file_path)
        if self.cancelled:
            self.complete.emit(False)
            return
        self.install_dependencies_and_cleanup()
        if not self.cancelled:
            self.complete.emit(True)

    def download_file(self,url,zip_file_path):
        """Downloads the file from the given URL."""
        try:
            with requests.get(url, stream=True) as response:
                response.raise_for_status()
                logger.debug(f"Started downloading from {url}.")
                self.process_download(response,zip_file_path)
        except requests.exceptions.RequestException as e:
            self.handle_network_error(e)

    def process_download(self, response,zip_file_path):
        """Processes the downloaded data."""

        with open(zip_file_path, 'wb') as f:
            self.write_chunks(response, f, zip_file_path)
        if not self.cancelled:
            self.status.emit("Downloading of zip file completed.")
            logger.info(f"Download completed for {self.url}.")


    def write_chunks(self, response, file_obj, zip_file_path):
        """Writes chunks of data to the file."""
        downloaded = 0
        self.last_logged_mb = 0

        for chunk in response.iter_content(chunk_size=4096):
            if self.cancelled:
                self.handle_download_cancellation(zip_file_path)
                return
            file_obj.write(chunk)
            downloaded += len(chunk)
            self.update_progress(downloaded)

    def update_progress(self, downloaded):
        """Updates the progress and emits status signals."""
        downloaded_mb = downloaded / (1024 * 1024)
        self.indeterminate.emit()
        if downloaded_mb - self.last_logged_mb >= 1:
            self.status.emit(f"Downloaded: {downloaded_mb:.2f} MB")
            logger.info(f"Downloaded: {downloaded_mb:.2f} MB")
            self.last_logged_mb = downloaded_mb

    def handle_download_cancellation(self, zip_file_path):
        """Handles download cancellation."""
        self.status.emit("Download canceled.")
        zip_file_path.unlink(missing_ok=True)
        logger.info(f"Download for {self.url} canceled by user.")

    def handle_network_error(self, error):
        """Handles network errors."""
        self.status.emit(f"Network error: {error}")
        logger.error(f"Network error while downloading {self.url}: {error}")

    def unzip_file(self):
        """Unzips the downloaded file and dynamically detects the extracted folder."""
        zip_file_name = self.url.split("/")[-1]
        zip_file_path = self.save_path / zip_file_name
        try:
            with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
                zip_ref.extractall(self.save_path)
                logger.info(f"File extracted to {self.save_path}.")
                self.status.emit("File unzipped successfully.")
            extracted_folder_path = self.detect_extracted_folder(zip_ref)
            logger.info(f"Detected extracted folder: {extracted_folder_path}")
        except (zipfile.BadZipFile, OSError) as e:
            if isinstance(e, zipfile.BadZipFile):
                error_message = "Failed to unzip: The file is corrupted or not a valid zip file."
            else:
                error_message = f"File system error while unzipping the file: {e}"
            self.status.emit(error_message)
            logger.error(f"Failed to unzip {zip_file_path}: {error_message}")
            raise

        return extracted_folder_path

    def detect_extracted_folder(self, zip_ref):
        """Detects the folder created during extraction."""
        folder_names = {Path(file.filename).parts[0] for file in zip_ref.infolist() if '/' in file.filename}

        if not folder_names:
            return self.save_path
        if len(folder_names) > 1:
            error_message = "Unexpected archive structure: multiple top-level folders detected."
            logger.error(error_message)
            raise ValueError(error_message)

        extracted_folder_name = folder_names.pop()
        extracted_folder_path = self.save_path / extracted_folder_name

        if not extracted_folder_path.is_dir():
            error_message = f"Expected extracted folder '{extracted_folder_name}' not found."
            logger.error(error_message)
            raise FileNotFoundError(error_message)

        return extracted_folder_path


    def install_dependencies_and_cleanup(self):
        """Installs dependencies and cleans up temporary files."""
        original_path = os.getcwd()
        zip_file_name = self.url.split("/")[-1]
        zip_file_path = self.save_path / zip_file_name
        self.status.emit("Unzipping downloaded file...")
        try:
            extracted_folder_path = self.unzip_file()
            if not extracted_folder_path:
                raise ValueError("Extraction failed. No extracted folder found.")

            os.chdir(extracted_folder_path)
            logger.info(f"Changed directory to {extracted_folder_path}")

            self.install_dependencies()
            logger.info("Dependencies installed successfully.")
        except Exception as e:
            error_message = f"Error during dependency installation: {e}"
            self.status.emit(error_message)
            logger.error(error_message)
        finally:
            os.chdir(original_path)
            logger.info(f"Returned to original directory: {original_path}")

            if zip_file_path.exists():
                zip_file_path.unlink(missing_ok=True)
                logger.info(f"Downloaded zip file {zip_file_path} has been deleted after extraction.")

    def cancel_download(self):
        self.cancelled = True
        logger.info(f"Cancellation requested for download of {self.url}.")


    def install_dependencies(self):
        """
        Downloads and installs Conda Forge on the system, and creates a Conda environment in the current directory.
        """
        if not self.cancelled:
            self.check_and_install_miniforge()
        if not self.cancelled:
            self.create_conda_environment()
            if not self.cancelled:
                logger.info("Environment created and dependencies installed successfully.")
                self.status.emit("Environment created and dependencies installed successfully.")

    def install_miniforge(self):
        """Download and install Miniforge if it is not already installed."""

        miniforge_url = "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
        installer_name = "Miniforge3-Linux-x86_64.sh"
        self.status.emit("Installing dependencies: downloading miniforge installer file...")
        logger.info(f"Downloading Miniforge installer...")
        installer_path = os.path.join(os.getcwd(), installer_name)
        self.download_file(miniforge_url,installer_path )
        logger.info(f"Miniforge installer downloaded: {installer_name}")

        logger.info("Running Miniforge installer...")
        subprocess.run(["bash", installer_path, "-b"], check=True)
        logger.info("Miniforge installed successfully.")
        self.status.emit("Installing dependencies: miniforge installed successfully.")

        logger.info("Removing Miniforge installer file...")
        os.remove(installer_path)
        logger.info(f"Miniforge installer file {installer_name} removed.")

    def check_and_install_miniforge(self):
        """Check if Miniforge (conda) is installed, and install it if necessary."""
        try:
            self.status.emit("Installing dependencies: checking if miniforge is installed...")
            logger.info("CHECKING FOR CONDA")
            subprocess.run(["which", "conda"], check=True, stdout=subprocess.PIPE)
            logger.info("Miniforge (conda) is already installed.")
        except subprocess.CalledProcessError:
            logger.warning("Miniforge (conda) is not installed. Installing Miniforge...")
            self.status.emit("Installing dependencies: installing miniforge...")
            self.install_miniforge()
        logger.info("CONDA FINISHED")

    def create_conda_environment(self):
        """
        Create a Conda environment in the current directory and log the process in real-time.
        :return:
        """

        def log_subprocess_output(pipe, level):
            """Log subprocess output in real-time."""
            for line in iter(pipe.readline, b''):
                logger.log(level, line.strip())
            pipe.close()

        self.status.emit("Installing dependencies: creating Conda environment...")
        env_dir = os.path.join(os.getcwd(), "snakemake")
        command = [
            "conda", "create", "--prefix", env_dir,
            "-c", "conda-forge", "-c", "bioconda",
            "minimap2", "snakemake", "last", "samtools", "-y"
        ]

        try:
            logger.info(f"Creating Conda environment at {env_dir}...")
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            log_subprocess_output(process.stdout, logging.INFO)
            log_subprocess_output(process.stderr, logging.ERROR)
            process.wait()
            if process.returncode == 0:
                self.status.emit("Installing dependencies: successfully created Conda environment.")
                logger.info("Conda environment created successfully.")
            else:
                raise subprocess.CalledProcessError(process.returncode, command)
        except subprocess.CalledProcessError as e:
            logger.error(f"Error occurred while creating environment: {e}")
            self.status.emit("Installing dependencies: error during creation of Conda environment.")

        self.status.emit("Installing dependencies: installing Python dependencies")
        pip_requirements_file = os.path.join(os.getcwd(), "requirements.txt")

        pip_command = [
            "conda", "run", "--prefix", env_dir,
            "python", "-m", "pip", "install",
            "--upgrade",
            "-r", pip_requirements_file
        ]

        if os.path.exists(pip_requirements_file):
            logger.info(f"Found requirements file: {pip_requirements_file}. Installing dependencies...")
            try:
                process = subprocess.Popen(pip_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                log_subprocess_output(process.stdout, logging.INFO)
                log_subprocess_output(process.stderr, logging.ERROR)
                process.wait()
                if process.returncode == 0:
                    logger.info("Python dependencies installed successfully from requirements.txt.")
                    self.status.emit(
                        "Installing dependencies: Python dependencies installed successfully from requirements.txt.")
                else:
                    raise subprocess.CalledProcessError(process.returncode, pip_command)
            except subprocess.CalledProcessError as e:
                logger.error(f"Error occurred while installing Python dependencies: {e}")
                self.status.emit("Installing dependencies: error during Python dependency installation.")
        else:
            logger.warning(
                f"No requirements.txt file found at {pip_requirements_file}. Skipping Python dependency installation.")

        list_packages_command = [
            "conda", "list", "--prefix", env_dir
        ]
        try:
            logger.info("Listing installed packages in the environment:")
            process = subprocess.Popen(list_packages_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            log_subprocess_output(process.stdout, logging.INFO)
            log_subprocess_output(process.stderr, logging.ERROR)
            process.wait()
            if process.returncode != 0:
                raise subprocess.CalledProcessError(process.returncode, list_packages_command)
        except subprocess.CalledProcessError as e:
            logger.error(f"Error occurred while listing installed packages: {e}")

class InstallerApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.installPathTextField = None
        self.back_button = None
        self.setFixedSize(700, 500)
        self.setWindowTitle("Oligo-bench Installer")

        self.download_path = None
        self.statusLabel = None
        self.is_downloading = None
        self.downloadButton = None
        self.progressBar = QProgressBar()
        self.installPathButton = None
        self.label = None
        self.layout = None
        self.downloadThread = None
        self.cancelButton = None
        self.isFinishedDownloading = False

        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)

        self.stackedWidget = QStackedWidget()
        self.initUI()
        self.apply_stylesheet()

    def apply_stylesheet(self):
        style_file = Path("styles/styles.qss")
        if style_file.exists():
            with open(style_file, 'r') as f:
                stylesheet = f.read()
                self.setStyleSheet(stylesheet)
        else:
            logger.warning(f"Style sheet {style_file} not found.")

    def initUI(self):
        main_layout = QVBoxLayout()

        # Navigation buttons
        nav_layout = QHBoxLayout()
        self.back_button = QPushButton("Back")
        self.back_button.clicked.connect(self.go_back)
        self.back_button.setEnabled(False)

        self.next_button = QPushButton("Next")
        self.next_button.clicked.connect(self.go_next)

        nav_layout.addWidget(self.back_button)
        nav_layout.addWidget(self.next_button)

        self.stackedWidget.addWidget(self.create_title_page())
        self.stackedWidget.addWidget(self.create_final_page())

        main_layout.addWidget(self.stackedWidget)
        main_layout.addLayout(nav_layout)
        self.central_widget.setLayout(main_layout)

    def create_title_page(self):
        page = QWidget()
        layout = QVBoxLayout()
        self.back_button.setVisible(False)
        title_label = QLabel("Welcome to Oligo-bench Installer")
        title_label.setFont(QFont("Arial", 16, QFont.Weight.Bold))
        title_label.setAlignment(Qt.AlignmentFlag.AlignCenter)

        layout.addWidget(title_label)
        page.setLayout(layout)
        return page

    def create_empty_page(self, page_number):
        page = QWidget()
        layout = QVBoxLayout()

        label = QLabel(f"Page {page_number}")
        label.setFont(QFont("Arial", 14))
        label.setAlignment(Qt.AlignmentFlag.AlignCenter)

        layout.addWidget(label)
        page.setLayout(layout)
        return page

    def create_final_page(self):
        page = QWidget()
        layout = QVBoxLayout()

        # Top Instruction Label
        self.topLabel = QLabel(
            "Setup will install Oligo-bench in the following folder. To install in a different folder, click Browse and select another folder. Click Next to continue.")
        self.topLabel.setWordWrap(True)
        layout.addWidget(self.topLabel)

        # Destination Folder Section
        pathLayout = QHBoxLayout()

        self.destLabel = QLabel("Destination Folder:")
        pathLayout.addWidget(self.destLabel, alignment=Qt.AlignmentFlag.AlignLeft)

        # Text field for installation path (initially set to current working directory)
        self.installPathTextField = QLineEdit("Choose directory...")  # Set the initial path to the current directory
        self.installPathTextField.setReadOnly(True)  # Make the text field read-only
        self.installPathTextField.setFixedWidth(300)
        pathLayout.addWidget(self.installPathTextField, alignment=Qt.AlignmentFlag.AlignCenter)

        # Browse button to open file explorer
        self.browseButton = QPushButton("Browse")
        self.browseButton.clicked.connect(self.choose_installation_path)
        pathLayout.addWidget(self.browseButton, alignment=Qt.AlignmentFlag.AlignRight)

        layout.addLayout(pathLayout)

        # Space Required and Available Section
        self.spaceRequiredLabel = QLabel(f"Space Required: 5 GB")
        layout.addWidget(self.spaceRequiredLabel)

        # Progress Bar (initially hidden)
        self.progressBar = QProgressBar()
        self.progressBar.setRange(0, 100)
        self.progressBar.setVisible(False)
        layout.addWidget(self.progressBar)

        # Status Label
        self.statusLabel = QLabel("")
        layout.addWidget(self.statusLabel)

        page.setLayout(layout)
        return page

    def go_next(self):
        current_index = self.stackedWidget.currentIndex()
        if current_index < self.stackedWidget.count() - 1:
            self.stackedWidget.setCurrentIndex(current_index + 1)
            self.back_button.setEnabled(True)
            self.back_button.setVisible(True)
        if current_index + 1 == self.stackedWidget.count() - 1:
            self.next_button.setText("Download")
            self.next_button.setEnabled(True)
        elif current_index + 1 == self.stackedWidget.count():
            if self.isFinishedDownloading:
                self.close()
            else:
                self.toggle_download()

    def go_back(self):
        current_index = self.stackedWidget.currentIndex()
        if self.is_downloading:
            QMessageBox.critical(self, "Error", "Downloading is in process...")
        else:
            if current_index > 0:
                self.stackedWidget.setCurrentIndex(current_index - 1)
                self.next_button.setEnabled(True)
                self.next_button.setText("Next")
            if current_index - 1 == 0:
                self.back_button.setEnabled(False)
                self.back_button.setVisible(False)

    def toggle_download(self):
        if self.is_downloading:
            self.cancel_download()
        else:
            self.start_download()

    def choose_installation_path(self):
        folder_dialog = QFileDialog(self)
        folder_dialog.setFileMode(QFileDialog.FileMode.Directory)
        folder_dialog.setOption(QFileDialog.Option.ShowDirsOnly)

        if folder_dialog.exec():
            self.download_path = folder_dialog.selectedFiles()[0]
            self.installPathTextField.setText(f"{self.download_path}")
            logger.info(f"Installation path chosen: {self.download_path}")

    def start_download(self):
        if not self.download_path:
            QMessageBox.critical(self, "Error", "Please choose an installation path first.")
            return

        self.progressBar.setVisible(True)

        self.downloadThread = DownloadThread("https://github.com/genomika-lt/oligo-bench/archive/refs/heads/main.zip", self.download_path)
        self.downloadThread.progress.connect(self.update_progress)
        self.downloadThread.status.connect(self.update_status)
        self.downloadThread.complete.connect(self.download_complete)
        self.downloadThread.indeterminate.connect(self.make_indeterminate)

        self.downloadThread.start()
        self.next_button.setText("Cancel Download")
        self.is_downloading = True

    def update_progress(self, progress):
        self.progressBar.setMaximum(100)
        self.progressBar.setValue(progress)
        logger.info(f"Progress bar updated to {progress}%")

    def make_indeterminate(self):
        self.progressBar.setMaximum(0)  # Indeterminate mode

    def update_status(self, status_message):
        self.statusLabel.setText(status_message)
        logger.info(status_message)

    def download_complete(self, success):
        self.progressBar.setMaximum(100)
        self.progressBar.setValue(100 if success else 0)
        if success:
            logger.info("Download process completed successfully.")
            self.update_status("Download process completed successfully.")
        else:
            logger.error("Download process failed.")
            self.update_status("Download process failed.")

        message = QMessageBox(self)
        message.setIcon(QMessageBox.Icon.Information if success else QMessageBox.Icon.Critical)
        message.setWindowTitle("Download Complete")
        message.setText("The download has completed successfully!" if success else "The download has failed.")
        message.setStandardButtons(QMessageBox.StandardButton.Ok)
        message.exec()

        self.next_button.setText("Finish")
        self.progressBar.setVisible(False)
        self.is_downloading = False
        self.isFinishedDownloading = True
        self.next_button.setEnabled(True)

    def cancel_download(self):
        if self.downloadThread:
            self.downloadThread.cancel_download()
            self.update_status("Download process was cancelled.")
            logger.info("Canceling download...")
            self.next_button.setText("Download")
            self.progressBar.setVisible(False)
            self.is_downloading = False

    def closeEvent(self, event):
        """
        Handles the window close event to manage unsaved
        changes before exiting the application.

        This method prompts the user to save any unsaved changes
        when the user attempts to close the window.
        It displays a confirmation dialog asking whether to save
        the changes, and if the user chooses to do so,
        it presents additional prompts to save the CSV and YAML
        files if applicable.

        :param event:
        :return:
        """
        # pylint: disable=C0103
        if self.isFinishedDownloading:
            event.accept()
        reply = QMessageBox.question(
            self,
            'Oligo-bench Setup',
            "Are you sure you want to quit Oligo-bench setup?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
            ,QMessageBox.StandardButton.No
        )
        if reply == QMessageBox.StandardButton.Yes:
            event.accept()
        elif reply == QMessageBox.StandardButton.No:
            event.ignore()

def main():
    app = QApplication(sys.argv)
    window = InstallerApp()
    window.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()