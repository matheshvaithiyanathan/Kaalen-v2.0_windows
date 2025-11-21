; -- Kaalen_Installer.iss --

[Setup]
AppName=Kaalen Data Viewer
AppVersion=1.0
AppPublisher=Mathesh Vaithiyanathan
DefaultDirName={autopf}\Kaalen
DefaultGroupName=Kaalen
AllowNoIcons=yes
OutputDir=Output_Installer
OutputBaseFilename=Kaalen_Setup_v1.0
Compression=lzma
SolidCompression=yes
SetupIconFile=icon.ico
WizardStyle=modern
RestartIfNeeded=yes
AppId=Kaalen_Data_Viewer
UninstallDisplayIcon={app}\Kaalen_App.exe
VersionInfoVersion=1.0.0
VersionInfoCompany=Mathesh Vaithiyanathan
PrivilegesRequired=admin

[Files]
Source: "dist\Kaalen_App\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\Kaalen Data Viewer"; Filename: "{app}\Kaalen_App.exe"; WorkingDir: "{app}"
Name: "{autodesktop}\Kaalen Data Viewer"; Filename: "{app}\Kaalen_App.exe"; WorkingDir: "{app}"; Tasks: desktopicon

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"

[Run]
Filename: "{app}\Kaalen_App.exe"; Description: "{cm:LaunchProgram,Kaalen Data Viewer}"; Flags: nowait postinstall skipifdoesntexist

[Registry]
Root: "HKLM"; Subkey: "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\Kaalen_Data_Viewer"; ValueType: string; ValueName: "DisplayName"; ValueData: "Kaalen Data Viewer"; Flags: uninsdeletekey
Root: "HKLM"; Subkey: "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\Kaalen_Data_Viewer"; ValueType: string; ValueName: "DisplayVersion"; ValueData: "1.0"
Root: "HKLM"; Subkey: "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\Kaalen_Data_Viewer"; ValueType: string; ValueName: "Publisher"; ValueData: "Mathesh Vaithiyanathan"
Root: "HKLM"; Subkey: "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\Kaalen_Data_Viewer"; ValueType: string; ValueName: "InstallLocation"; ValueData: "{app}"
Root: "HKLM"; Subkey: "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\Kaalen_Data_Viewer"; ValueType: string; ValueName: "UninstallString"; ValueData: "{uninstallexe}"
Root: "HKLM"; Subkey: "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\Kaalen_Data_Viewer"; ValueType: string; ValueName: "DisplayIcon"; ValueData: "{app}\Kaalen_App.exe"
Root: "HKLM"; Subkey: "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\Kaalen_Data_Viewer"; ValueType: dword; ValueName: "NoModify"; ValueData: 1
Root: "HKLM"; Subkey: "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\Kaalen_Data_Viewer"; ValueType: dword; ValueName: "NoRepair"; ValueData: 1

[UninstallDelete]
Type: filesandordirs; Name: "{app}\*"
