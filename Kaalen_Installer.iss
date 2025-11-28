; -- Kaalen_Installer_v2.0.iss --
; The issue was caused by the manual [Registry] section creating a second, conflicting entry.

[Setup]
; --- Application Identity ---
AppName=Kaalen Version 2.0  ; The single, clean name for the Add/Remove Programs list
AppVersion=2.0
AppPublisher=Mathesh Vaithiyanathan
VersionInfoCompany=Mathesh Vaithiyanathan (Author)

; --- CRITICAL FIX: The AppId ensures only ONE uninstall entry is created. ---
; *** Do NOT change this GUID in future versions (e.g., v2.1, v3.0). ***
AppId={{50E0D2F8-3A7B-46C9-A1C8-710E1C92E152} 

; --- Installation Paths ---
DefaultDirName={autopf}\Kaalen
DefaultGroupName=Kaalen
AllowNoIcons=yes
SetupIconFile=icon.ico
WizardStyle=modern

; --- Output Settings ---
OutputDir=Output_Installer
OutputBaseFilename=Kaalen_Setup_v2.0
Compression=lzma
SolidCompression=yes

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

; --- IMPORTANT: The entire [Registry] section has been REMOVED ---
; Inno Setup automatically handles the creation of the uninstall registry key 
; when the AppId is defined, solving your issue.
