; -- Kaalen_Installer.iss --

[Setup]
; IMPORTANT: Change the AppName, AppVersion, AppPublisher, and icon path as needed
AppName=Kaalen 
AppVersion=2.0
AppPublisher=Mathesh Vaithiyanathan
VersionInfoCompany=Mathesh Vaithiyanathan (Author)
DefaultDirName={autopf}\Kaalen
DefaultGroupName=Kaalen
AllowNoIcons=yes
OutputDir=Output_Installer
OutputBaseFilename=Kaalen_Setup_v2.0
Compression=lzma
SolidCompression=yes
SetupIconFile=icon.ico 
WizardStyle=modern
; Removed invalid directive: RestartIfNeeded=yes 
; If you need restart functionality, use one of these instead:
; RestartApplications=yes
; Or for more control: AlwaysRestart=no

[Files]
; This command copies all files (recursively) from the PyInstaller build directory
Source: "dist\Kaalen_App\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
; Creates a Start Menu shortcut
Name: "{group}\Kaalen"; Filename: "{app}\Kaalen_App.exe"; IconFilename: "{app}\Kaalen_App.exe"; WorkingDir: "{app}"
; Creates an optional Desktop shortcut
Name: "{autodesktop}\Kaalen"; Filename: "{app}\Kaalen_App.exe"; IconFilename: "{app}\Kaalen_App.exe"; WorkingDir: "{app}"; Tasks: desktopicon

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}";

[Run]
Filename: "{app}\Kaalen_App.exe"; Description: "{cm:LaunchProgram,Kaalen}"; Flags: nowait postinstall skipifdoesntexist

[Registry]
; Adds the application to Windows' Add/Remove Programs list
Root: "HKLM"; Subkey: "Software\Microsoft\Windows\CurrentVersion\Uninstall\Kaalen_App"; ValueType: "string"; ValueName: "DisplayName"; ValueData: "Kaalen Data Viewer";
Root: "HKLM"; Subkey: "Software\Microsoft\Windows\CurrentVersion\Uninstall\Kaalen_App"; ValueType: "string"; ValueName: "DisplayVersion"; ValueData: "1.0";
Root: "HKLM"; Subkey: "Software\Microsoft\Windows\CurrentVersion\Uninstall\Kaalen_App"; ValueType: "string"; ValueName: "Publisher"; ValueData: "Mathesh Vaithiyanathan";
Root: "HKLM"; Subkey: "Software\Microsoft\Windows\CurrentVersion\Uninstall\Kaalen_App"; ValueType: "string"; ValueName: "InstallLocation"; ValueData: "{app}";
Root: "HKLM"; Subkey: "Software\Microsoft\Windows\CurrentVersion\Uninstall\Kaalen_App"; ValueType: "string"; ValueName: "UninstallString"; ValueData: "{uninstallexe}";
