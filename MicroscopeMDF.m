% Most Software Machine Data File

%% scanimage.SI (ScanImage)

% Global microscope properties
objectiveResolution = 10;     % Resolution of the objective in microns/degree of scan angle

% Data file location

% Custom Scripts
startUpScript = '';     % Name of script that is executed in workspace 'base' after scanimage initializes
shutDownScript = '';     % Name of script that is executed in workspace 'base' after scanimage exits

fieldCurvatureZs = [];     % Field curvature for mesoscope
fieldCurvatureRxs = [];     % Field curvature for mesoscope
fieldCurvatureRys = [];     % Field curvature for mesoscope
fieldCurvatureTip = 0;     % Field tip for mesoscope
fieldCurvatureTilt = 0;     % Field tilt for mesoscope
useJsonHeaderFormat = false;     % Use JSON format for TIFF file header

%% scanimage.components.CoordinateSystems (SI CoordinateSystems)
% SI Coordinate System Component.
classDataFileName = 'default-CoordinateSystems_classData.mat';     % File containing the previously generated alignment data corresponding to the currently installed objective, SLM, scanners, etc.

%% scanimage.components.Motors (SI Motors)
% SI Stage/Motor Component.
motorXYZ = {'' '' ''};     % Defines the motor for ScanImage axes X Y Z.
motorAxisXYZ = [1 2 3];     % Defines the motor axis used for Scanimage axes X Y Z.
scaleXYZ = [1 1 1];     % Defines scaling factors for axes.
backlashCompensation = [0 0 0];     % Backlash compensation in um (positive or negative)
moveTimeout_s = 10;     % Move timeout in seconds

%% scanimage.components.Photostim (SI Photostim)
photostimScannerName = '';     % Name of scanner (from first MDF section) to use for photostimulation. Must be a linear scanner
stimTriggerTerm = 1;     % Specifies the channel that should be used to trigger a stimulation. This a triggering port name such as D2.1 for vDAQ or PFI1 for the auxiliary IO board of an NI LinScan system.

% Monitoring DAQ AI channels
BeamAiId = [];     % AI channel to be used for monitoring the Pockels cell output

loggingStartTrigger = '';     % PFI line to which start trigger for logging is wired to photostim board. Leave empty for automatic routing via PXI bus

stimActiveOutputChannel = '';     % Digital terminal on stim board to output stim active signal. (e.g. on vDAQ: 'D2.6' on NI-DAQ hardware: '/port0/line0'
beamActiveOutputChannel = '';     % Digital terminal on stim board to output beam active signal. (e.g. on vDAQ: 'D2.7' on NI-DAQ hardware: '/port0/line1'
slmTriggerOutputChannel = '';     % Digital terminal on stim board to trigger SLM frame flip. (e.g. on vDAQ: 'D2.5' on NI-DAQ hardware: '/port0/line2'

%% scanimage.components.scan2d.LinScan (gals)
deviceNameAcq = 'Device1';     % string identifying NI DAQ board for PMT channels input
deviceNameAux = 'Device2';     % string identifying NI DAQ board for outputting clocks. leave empty if unused. Must be a X-series board

externalSampleClock = false;     % Logical: use external sample clock connected to the CLK IN terminal of the FlexRIO digitizer module
externalSampleClockRate = 8e+07;     % [Hz]: nominal frequency of the external sample clock connected to the CLK IN terminal (e.g. 80e6); actual rate is measured on FPGA

% Optional
channelsInvert = [false false false false];     % scalar or vector identifiying channels to invert. if scalar, the value is applied to all channels

xGalvo = 'Galvo X';     % x Galvo device name
yGalvo = 'Galvo Y';     % y Galvo device name
fastZs = {};     % fastZ device names
beams = {};     % fastZ device names
shutters = {'Laser shutter'};     % shutter device names

referenceClockIn = '';     % one of {'',PFI14} to which 10MHz reference clock is connected on Aux board. Leave empty for automatic routing via PXI bus
enableRefClkOutput = false;     % Enables/disables the export of the 10MHz reference clock on PFI14

% Acquisition
channelIDs = [0 1 2 3];     % Array of numeric channel IDs for PMT inputs. Leave empty for default channels (AI0...AIN-1)

% Advanced/Optional:
stripingEnable = true;     % enables/disables striping display
stripingMaxRate = 10;     % [Hz] determines the maximum display update rate for striping
maxDisplayRate = 30;     % [Hz] limits the maximum display rate (affects frame batching)
internalRefClockSrc = '';     % Reference clock to use internally
internalRefClockRate = [];     % Rate of reference clock to use internally
secondaryFpgaFifo = false;     % specifies if the secondary fpga fifo should be used

% Laser Trigger
LaserTriggerPort = '';     % Port on FlexRIO AM digital breakout (DIO0.[0:3]) or digital IO DAQ (PFI[0:23]) where laser trigger is connected.
LaserTriggerFilterTicks = 0;
LaserTriggerSampleMaskEnable = false;
LaserTriggerSampleWindow = [0 1];

% Calibration data
scannerToRefTransform = [1 0 0;0 1 0;0 0 1];

%% dabs.generic.GalvoPureAnalog (Galvo Y)
AOControl = '/Device1/AO1';     % control terminal  e.g. '/vDAQ0/AO0'
AOOffset = '';     % control terminal  e.g. '/vDAQ0/AO0'
AIFeedback = '';     % feedback terminal e.g. '/vDAQ0/AI0'

angularRange = 40;     % total angular range in optical degrees (e.g. for a galvo with -20..+20 optical degrees, enter 40)
voltsPerOpticalDegrees = 0.5;     % volts per optical degrees for the control signal
voltsOffset = 0;     % voltage to be added to the output
parkPosition = 20;     % park position in optical degrees
slewRateLimit = Inf;     % Slew rate limit of the analog output in Volts per second

% Calibration settings
feedbackVoltLUT = zeros(0,2);     % [Nx2] lut translating feedback volts into position volts
offsetVoltScaling = 1;     % scalar factor for offset volts

%% dabs.generic.GalvoPureAnalog (Galvo X)
AOControl = '/Device1/AO0';     % control terminal  e.g. '/vDAQ0/AO0'
AOOffset = '';     % control terminal  e.g. '/vDAQ0/AO0'
AIFeedback = '';     % feedback terminal e.g. '/vDAQ0/AI0'

angularRange = 40;     % total angular range in optical degrees (e.g. for a galvo with -20..+20 optical degrees, enter 40)
voltsPerOpticalDegrees = 0.5;     % volts per optical degrees for the control signal
voltsOffset = 0;     % voltage to be added to the output
parkPosition = 20;     % park position in optical degrees
slewRateLimit = Inf;     % Slew rate limit of the analog output in Volts per second

% Calibration settings
feedbackVoltLUT = zeros(0,2);     % [Nx2] lut translating feedback volts into position volts
offsetVoltScaling = 1;     % scalar factor for offset volts

%% dabs.generic.DigitalShutter (Laser shutter)
DOControl = '/Device2/port0/line7';     % control terminal  e.g. '/vDAQ0/DIO0'
invertOutput = false;     % invert output drive signal to shutter
openTime_s = 0.5;     % settling time for shutter in seconds
shutterTarget = 'Excitation';     % one of {', 'Excitation', 'Detection'}

%% dabs.generic.PMTAnalog (SHG PMT)
AOGain = '/Device2/AO0';     % <optional> resource name of the analog output channel that controls the PMT gain (e.g. /vDAQ0/AO0)
AOOutputRange = [0 10];     % <required if AOGain is defined> array of 1x2 numeric array specifying the minimum and maximum analog output voltage on the DAQ board that controls the PMT gain.
SupplyVoltageRange = [0 1250];     % <required if AOGain is defined> array of 1x2 specifying the minimum and maximum for the PMT power supply in Volts.

DOPower = '';     % <optional> resource name of the digital output channel that switches the PMT on/off (e.g. /vDAQ0/D0.0)
DITripDetect = '';     % <optional> resource name of the analog output channel that controls the PMT gain (e.g. /vDAQ0/D0.1)
DOTripReset = '';     % <optional> resource name of the analog output channel that controls the PMT gain (e.g. /vDAQ0/D0.2)

% Calibration settings
wavelength_nm = 509;     % wavelength in nanometer
autoOn = true;     % powers the PMT automatically on for the duration of a scan
gain_V = 600;     % PMT power supply voltage

%% dabs.legacy.motor.LegacyMotor (Prior stage)
% Motor used for X/Y/Z motion, including stacks.

controllerType = 'Prior ProScan III';     % If supplied, one of {'sutter.mp285', 'sutter.mpc200', 'thorlabs.mcm3000', 'thorlabs.mcm5000', 'scientifica', 'pi.e665', 'pi.e816', 'npoint.lc40x', 'bruker.MAMC'}.
comPort = 4;     % Integer identifying COM port for controller, if using serial communication
customArgs = {};     % Additional arguments to stage controller. Some controller require a valid stageType be specified
invertDim = '';     % string with one character for each dimension specifying if the dimension should be inverted. '+' for normal, '-' for inverted
positionDeviceUnits = [];     % 1xN array specifying, in meters, raw units in which motor controller reports position. If unspecified, default positionDeviceUnits for stage/controller type presumed.
velocitySlow = [];     % Velocity to use for moves smaller than motorFastMotionThreshold value. If unspecified, default value used for controller. Specified in units appropriate to controller type.
velocityFast = [];     % Velocity to use for moves larger than motorFastMotionThreshold value. If unspecified, default value used for controller. Specified in units appropriate to controller type.
moveCompleteDelay = [];     % Delay from when stage controller reports move is complete until move is actually considered complete. Allows settling time for motor
moveTimeout = '';     % Default: 2s. Fixed time to wait for motor to complete movement before throwing a timeout error
moveTimeoutFactor = '';     % (s/um) Time to add to timeout duration based on distance of motor move command

%% dabs.generic.MotorizedHalfWavePlate (HWP power)
rotationStage = 'HWP K Cube motor';     % Motor Controller user-assigned device name  e.g. 'Rotation Stage'
motorAxis = 1;     % Number of the axis for on the motor controller for controlling the Rotation Stage
AIFeedback = '/Device2/AI3';     % feedback terminal e.g. '/vDAQ0/AI0'

moveTimeout_s = 100;     % move timeout in seconds
outputRange_deg = [0 90];     % Control angular range in degrees
devUnitsPerDegree = 1;     % Ratio of stage units per degree of rotation
feedbackUsesRejectedLight = false;     % Indicates if photodiode is in rejected path of beams modulator.
calibrationOpenShutters = {'Laser shutter'};     % List of shutters to open during the calibration. (e.g. {'Shutter1' 'Shutter2'}

powerFractionLimit = 1;     % Maximum allowed power fraction (between 0 and 1)

% Calibration data
powerFraction2ModulationAngleLut = [0 15;0.0292398 18.75;0.0890953 22.5;0.179051 26.25;0.291108 30;0.416323 33.75;0.54859 37.5;0.676385 41.25;0.791452 45;0.889061 48.75;0.95872 52.5;0.994582 56.25;1 60];
powerFraction2PowerWattLut = [0 0;1 0.15];
powerFraction2FeedbackVoltLut = [0 0.00477578;1 0.753378];
feedbackOffset_V = 0;

% Calibration settings
calibrationNumPoints = 25;     % number of equidistant points to measure within the angular output range
calibrationAverageSamples = 5;     % per calibration point, average N analog input samples. This helps to reduce noise
calibrationNumRepeats = 1;     % number of times to repeat the calibration routine. the end result is the average of all calibration runs
calibrationMotorSettlingTime_s = 5;     % pause between measurement points. this allows the beam modulation to settle

%% dabs.thorlabs.KinesisMotor (HWP K Cube motor)
serial = '27254531';     % Serial of the thorlabs stage e.g. '45155204'
kinesisInstallDir = 'C:\Program Files\Thorlabs\Kinesis';     % Path to Thorlabs Kinesis installation
homingTimeout_s = 20;     % Timeout for homing move in seconds
units = 'deg';     % Units for the device. One of {'um' 'deg'}
startupSettingsMode = 'UseDeviceSettings';     % Settings the device should use. One of {"Reset RealToDeviceUnit","Use Device Settings","Use File Settings","Use Configured Settings"}

%% dabs.generic.PMTAnalog (TPEF 1 PMT)
AOGain = '/Device2/AO1';     % <optional> resource name of the analog output channel that controls the PMT gain (e.g. /vDAQ0/AO0)
AOOutputRange = [0 10];     % <required if AOGain is defined> array of 1x2 numeric array specifying the minimum and maximum analog output voltage on the DAQ board that controls the PMT gain.
SupplyVoltageRange = [0 1250];     % <required if AOGain is defined> array of 1x2 specifying the minimum and maximum for the PMT power supply in Volts.

DOPower = '';     % <optional> resource name of the digital output channel that switches the PMT on/off (e.g. /vDAQ0/D0.0)
DITripDetect = '';     % <optional> resource name of the analog output channel that controls the PMT gain (e.g. /vDAQ0/D0.1)
DOTripReset = '';     % <optional> resource name of the analog output channel that controls the PMT gain (e.g. /vDAQ0/D0.2)

% Calibration settings
wavelength_nm = 509;     % wavelength in nanometer
autoOn = true;     % powers the PMT automatically on for the duration of a scan
gain_V = 600;     % PMT power supply voltage

%% dabs.generic.PMTAnalog (TPEF 2 PMT)
AOGain = '/Device2/AO2';     % <optional> resource name of the analog output channel that controls the PMT gain (e.g. /vDAQ0/AO0)
AOOutputRange = [0 5];     % <required if AOGain is defined> array of 1x2 numeric array specifying the minimum and maximum analog output voltage on the DAQ board that controls the PMT gain.
SupplyVoltageRange = [0 1250];     % <required if AOGain is defined> array of 1x2 specifying the minimum and maximum for the PMT power supply in Volts.

DOPower = '';     % <optional> resource name of the digital output channel that switches the PMT on/off (e.g. /vDAQ0/D0.0)
DITripDetect = '';     % <optional> resource name of the analog output channel that controls the PMT gain (e.g. /vDAQ0/D0.1)
DOTripReset = '';     % <optional> resource name of the analog output channel that controls the PMT gain (e.g. /vDAQ0/D0.2)

% Calibration settings
wavelength_nm = 509;     % wavelength in nanometer
autoOn = true;     % powers the PMT automatically on for the duration of a scan
gain_V = 600;     % PMT power supply voltage

