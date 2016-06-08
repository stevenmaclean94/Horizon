function varargout = player(varargin)
% PLAYER M-file for player.fig
%      PLAYER, by itself, creates a new PLAYER or raises the existing
%      singleton*.
%
%      H = PLAYER returns the handle to a new PLAYER or the handle to
%      the existing singleton*.
%
%      PLAYER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLAYER.M with the given input arguments.
%
%      PLAYER('Property','Value',...) creates a new PLAYER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before player_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to player_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help player

% Last Modified by GUIDE v2.5 14-Aug-2006 21:41:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @player_OpeningFcn, ...
                   'gui_OutputFcn',  @player_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before player is made visible.
function player_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to player (see VARARGIN)

% Choose default command line output for player
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

handles.fig = figure;
set(handles.fig, 'renderer', 'opengl');

hold on
grid on

[playIcon, map] = imread('play.bmp');
playIcon = ind2rgb(playIcon, map);
set(handles.play, 'cdata', playIcon)
% Load Scedule Data
handles.schedData = xlsread('schedData');
% Plot position of Target
plot3(handles.schedData(1,6), handles.schedData(1,7), handles.schedData(1,8), 'r*')
plot3(handles.schedData(1,18), handles.schedData(1,19), handles.schedData(1,20), 'g*')

% Load s/c position and Velocity Data
handles.propData = xlsread('propData');
% Plot s/c trajectory
plot3(handles.propData(:,2), handles.propData(:,3), handles.propData(:,4))

% Define Constants
handles.Re = 6378;
handles.mu = 398601;

% Create Earth
[sx sy sz] = sphere(360);
sx = sx*handles.Re;
sy = sy*handles.Re;
sz = sz*handles.Re;
I = imread('earthbare.jpg');
handles.earth = surf(sx,sy,sz);
% Rotate Earth so aligned with GMT
rotate(handles.earth, [0 0 1], handles.schedData(2,1));
set(handles.earth, 'edgeColor', 'none', 'CData', I, 'Facecolor', 'texturemap');

% Set lighting
handles.lig = camlight('infinite');
set(handles.lig, 'position', [0 1 0]);
lighting gouraud

xlabel 'x'
ylabel 'y'
zlabel 'z'
cameratoolbar('show')
cameratoolbar('ResetCameraAndSceneLight')

r = handles.propData(1,2:4)';
v = handles.propData(1,5:7)';

% Create Horizon View Cone
r_n = norm(r);
betaH = asin(handles.Re/r_n);
Dh = r_n*cos(betaH);
y = Dh * cos(betaH);
x = Dh * sin(betaH);
[cx cy cz] = cylinder(0:x, 36);
cz = cz * y;
handles.viewCone = surf(cx,cy,cz, 'edgecolor', 'none', 'facecolor', [.6 .1 .1], 'facealpha', .4);
% rotate cone and move into proper location
rotVec = cross(r/r_n, [0 0 1]);
rotAngle = 180 - 180/pi*acos(r'*[0; 0; 1]/r_n);
rotate(handles.viewCone, rotVec, rotAngle)
translate(handles.viewCone, r(1), r(2), r(3));

% create s/c
handles.aeolis = aeolis;
translate(handles.aeolis, r(1), r(2), r(3));
rotate(handles.aeolis, v, 90, r);

% calculate angular momentum vector
h = cross(r, v);
handles.h_hat = h/norm(h);

% calculate the orbital period
handles.T = 2*pi/sqrt(handles.mu)*(handles.Re+300)^(3/2);

r = interp1(handles.propData(:,1), handles.propData(:,2:4), 0)';
handles.p = r*1.01 + 20*handles.h_hat;
handles.u = r;
alpha = asin(handles.Re/norm(r)) + 2*pi/180;
norm_l = sqrt(norm(r)^2 - handles.Re^2);
lr = norm_l*cos(alpha);
lv = norm_l*sin(alpha);
v = interp1(handles.propData(:,1), handles.propData(:,5:7), 0)';
%propData(i,5:7)';
l = lr*(-r/norm(r)) + lv*(v/norm(v));
% camera target location
handles.targ = r + l;
u = r;
camup(handles.u);
campos(handles.p);
camtarget(handles.targ)

evalin('base', ['simStart = ',num2str(handles.propData(1,1)),';']);
evalin('base', ['simEnd = ',num2str(handles.propData(end,1)),';']);
evalin('base', ['currentSimTime = 0;']);
evalin('base', 'simStep = 100;');
evalin('base', 'simOn = 0;');

set(handles.slider, 'min', evalin('base', 'simStart;'));
set(handles.slider, 'max', evalin('base', 'simEnd;'));
set(handles.time, 'string', '0');

% UIWAIT makes player wait for user response (see UIRESUME)
% uiwait(handles.player_figure);

guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = player_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider_Callback(hObject, eventdata, handles)
% hObject    handle to slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

figure(handles.fig);

newSimTime = get(hObject, 'value');
paint_sim(handles, newSimTime, evalin('base', 'currentSimTime;'));
evalin('base', ['currentSimTime = ', num2str(newSimTime),';']);
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in rewind.
function rewind_Callback(hObject, eventdata, handles)
% hObject    handle to rewind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.simStep = handles.simStep*2;

guidata(hObject, handles)

% --- Executes on button press in frame_reverse.
function frame_reverse_Callback(hObject, eventdata, handles)
% hObject    handle to frame_reverse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure(handles.fig)

% If playing, stop!
if evalin('base', 'simOn') == 1
    play_Callback(handles.play, eventdata, handles)
end

% Advance ONE frame (simStep)
simTime = evalin('base', 'currentSimTime - simStep;');
if (simTime > evalin('base', 'simEnd;'))
    simTime = evalin('base', 'simEnd;');
end
paint_sim(handles, simTime, simTime + evalin('base', 'simStep;'));
evalin('base', ['currentSimTime = ', num2str(simTime),';']);
guidata(hObject, handles)

% --- Executes on button press in fast_forward.
function fast_forward_Callback(hObject, eventdata, handles)
% hObject    handle to fast_forward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.simStep = handles.simStep*2;
guidata(hObject, handles)
play_Callback(handles.play, eventdata, handles);

% --- Executes on button press in frame_advance.
function frame_advance_Callback(hObject, eventdata, handles)
% hObject    handle to frame_advance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure(handles.fig)

% If playing, stop!
if evalin('base', 'simOn') == 1
    play_Callback(handles.play, eventdata, handles)
end

% Advance ONE frame (simStep)
simTime = evalin('base', 'currentSimTime + simStep;');
if (simTime > evalin('base', 'simEnd;'))
    simTime = evalin('base', 'simEnd;');
end
paint_sim(handles, simTime, simTime - evalin('base', 'simStep;'));
evalin('base', ['currentSimTime = ', num2str(simTime),';']);
guidata(hObject, handles)


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes when user attempts to close player_figure.
function player_figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to player_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(handles.fig);
% Hint: delete(hObject) closes the figure
delete(hObject);

function paint_sim(handles, newSimTime, oldSimTime)

% determine how to translate s/c
dt = newSimTime - oldSimTime;%(t - t_old);

r_old = interp1(handles.propData(:,1), handles.propData(:,2:4), oldSimTime,'spline')';
%propData(i-dir,2), propData(i-dir,3), propData(i-dir,4)]';
r = interp1(handles.propData(:,1), handles.propData(:,2:4), newSimTime,'spline')';
%[propData(i,2), propData(i,3), propData(i,4)]';
d_r_sat = (r - r_old);
% determine s/c and view cone rotation angle
%    t_old = propData(i - dir,1);
%    t = propData(i,1);
ang = 360*dt/handles.T;

% Move s/c and view cone
translate(handles.aeolis, d_r_sat(1), d_r_sat(2), d_r_sat(3));
%     translate(panelL, d_r_sat(1), d_r_sat(2), d_r_sat(3));
%     translate(panelR, d_r_sat(1), d_r_sat(2), d_r_sat(3));
translate(handles.viewCone, d_r_sat(1), d_r_sat(2), d_r_sat(3));
% determine horizon angle for viewing
% alpha = asin(handles.Re/norm(r)) + 2*pi/180;
% norm_l = sqrt(norm(r)^2 - handles.Re^2);
% lr = norm_l*cos(alpha);
% lv = norm_l*sin(alpha);
% v = interp1(handles.propData(:,1), handles.propData(:,5:7), newSimTime)';
% %propData(i,5:7)';
% l = lr*(-r/norm(r)) + lv*(v/norm(v));
% % camea target location
% targ = r + l;

% rotate s/c and view cone
rotate(handles.aeolis, handles.h_hat, ang, r);
rotate(handles.viewCone, handles.h_hat, ang, r);

% camera position vector
%p = r*1.01 + 20*handles.h_hat;
targ = camtarget';
t_h = line([0 targ(1)], [0 targ(2)], [0 targ(3)], 'visible', 'off');
rotate(t_h, handles.h_hat, ang);
t_x = get(t_h, 'xdata');
t_y = get(t_h, 'ydata');
t_z = get(t_h, 'zdata');
targ = [t_x(2), t_y(2), t_z(2)];

u = camup';
u_h = line([0 u(1)], [0 u(2)], [0 u(3)], 'visible', 'off');
rotate(u_h, handles.h_hat, ang);
u_x = get(u_h, 'xdata');
u_y = get(u_h, 'ydata');
u_z = get(u_h, 'zdata');
u = [u_x(2), u_y(2), u_z(2)];

p = campos';
p_h = line([0 p(1)], [0 p(2)], [0 p(3)], 'visible', 'off');
rotate(p_h, handles.h_hat, ang);
p_x = get(p_h, 'xdata');
p_y = get(p_h, 'ydata');
p_z = get(p_h, 'zdata');
p = [p_x(2), p_y(2), p_z(2)];

campos(p);
camup(u);
%camva(30);
camtarget(targ);

evalin('base', ['currentSimTime = ', num2str(newSimTime),';']);
set(handles.time, 'string', newSimTime);
set(handles.slider, 'value', newSimTime);
guidata(handles.fig, handles);

drawnow


% --- Executes on button press in play.
function play_Callback(hObject, eventdata, handles)
% hObject    handle to play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure(handles.fig)

if evalin('base','simOn') == 0
    evalin('base','simOn = 1');
elseif evalin('base','simOn') == 1
    evalin('base', 'simOn = 0');
end

% Play Simulation
set(hObject, 'string', '"')
cst = evalin('base', 'currentSimTime;');
ss = evalin('base', 'simStep;');
se = evalin('base', 'simEnd;');
for simTime = cst+ss:ss:se
    if evalin('base', 'simOn') == 1
        paint_sim(handles, simTime, simTime - ss);
    else
        break;
    end
end

evalin('base', ['currentSimTime = ', num2str(simTime), '- simStep;']);
guidata(hObject, handles);
set(hObject, 'string', '>')