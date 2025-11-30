; -- Kaalen_Installer.iss --
; Fixes:
; 1. Uses a static AppId to prevent duplicate "Programs and Features" entries.
; 2. Uses [InstallDelete] to forcefully remove the old, broken registry entry ("Kaalen Data Viewer 1.0").

[Setup]
; --- Application Identity ---
AppName=Kaalen version 2.0         ; The desired display name for the Control Panel
AppVersion=2.0
AppPublisher=Mathesh Vaithiyanathan
VersionInfoCompany=Mathesh Vaithiyanathan (Author)

; *** CRITICAL FIX: The AppId must be static across all versions (v1.0, v2.0, etc.) ***
; This ensures that upgrading replaces the previous entry instead of adding a new one.
AppId={{50E0D2F8-3A7B-46C9-A1C8-710E1C92E152} 

; --- Installation Paths ---
DefaultDirName={autopf}\Kaalen
DefaultGroupName=Kaalen
AllowNoIcons=yes
OutputDir=Output_Installer
OutputBaseFilename=Kaalen_Setup_v2.0
Compression=lzma
SolidCompression=yes
SetupIconFile=icon.ico
WizardStyle=modern

; --- CRITICAL FIX: CLEANUP FOR BROKEN V1.0 REGISTRY ENTRY ---
; This targets the old, manual registry subkey that was creating the leftover entry 
; in Programs and Features (Kaalen Data Viewer 1.0) and deletes it during installation.
[InstallDelete]
; Target HKLM (Local Machine) where your previous script placed the entry
Type: regkey; Name: "HKEYLM\SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\Kaalen_App"

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

; *** IMPORTANT ***
; The manual [Registry] section from your previous script that created the duplicate
; "Kaalen Data Viewer 1.0" entry is now completely removed, as AppId handles this automatically.
