function varargout = ArcasGui(varargin)
% ARCASGUI MATLAB code for ArcasGui.fig
%      ARCASGUI, by itself, creates a new ARCASGUI or raises the existing
%      singleton*.
%
%      H = ARCASGUI returns the handle to a new ARCASGUI or the handle to
%      the existing singleton*.
%
%      ARCASGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ARCASGUI.M with the given input arguments.
%
%      ARCASGUI('Property','Value',...) creates a new ARCASGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ArcasGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ArcasGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ArcasGui

% Last Modified by GUIDE v2.5 26-May-2014 10:24:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ArcasGui_OpeningFcn, ...
    'gui_OutputFcn',  @ArcasGui_OutputFcn, ...
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


% --- Executes just before ArcasGui is made visible.
function ArcasGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ArcasGui (see VARARGIN)

% Choose default command line output for ArcasGui
handles.output = hObject;

% Disable analyzius buttons
set(handles.singleTrial_pushbutton,'Enable','off');
set(handles.runSimulation_pushbutton,'Enable','off');
set(handles.showHistogram_pushbutton,'Enable','off');

% Disable controls for thresholding
setThresholdButtons(handles.redThreshold_checkbox, [],handles);
setThresholdButtons(handles.blueThreshold_checkbox, [],handles);
setThresholdButtons(handles.greenThreshold_checkbox, [],handles);

% Disable thresholding checkboxes
set(handles.redThreshold_checkbox,'Enable','inactive');
set(handles.greenThreshold_checkbox,'Enable','inactive');
set(handles.blueThreshold_checkbox,'Enable','inactive');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ArcasGui wait for user response (see UIRESUME)
% uiwait(handles.ArcasGui_figure);


% --- Outputs from this function are returned to the command line.
function varargout = ArcasGui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadImage_pushbutton.
function loadImage_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to loadImage_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Prompt user to select image, try to import selected image, if fail prompt
% user again (or user can hit cancel to exit loop
DONE=0;
while ~DONE
    [FileName,PathName] = uigetfile('*.*','Select an Image for Simulation');
    if PathName==0; return; end
    
    try
        img = imread(strcat(PathName,FileName)); DONE = 1;
    catch ME; h = errordlg('File cannot be read by MATLAB as Image');
        uiwait(h);
    end
end
% Store path data for image
setappdata(handles.ArcasGui_figure, 'img_path', strcat(PathName,FileName));



% Display image
imshow(img, 'Parent', handles.img_axes);

elem = @(x) x(:);
% Actions taken depend on nature on input image
if islogical(img) && size(img,3)==1
    % Assign image to bw image
    setappdata(handles.ArcasGui_figure, 'bw_img',img);
    % Convert image to RGB
    setappdata(handles.ArcasGui_figure, 'img',cat(3,im2uint8(img),im2uint8(img),...
        im2uint8(img)));
    % Enable trial buttons
    trial_params_changed(handles);
% Grayscale image
elseif size(img,3)==1
    % Assign image to bw image
    setappdata(handles.ArcasGui_figure, 'bw_img',[]);
    % Convert image to RGB
    setappdata(handles.ArcasGui_figure, 'img',cat(3,im2uint8(img),im2uint8(img),...
        im2uint8(img)));
    % Enable trial buttons
    disable_trial_buttons(handles);   
% If all channels of image equal, 2 unique values, then its a 3 channel bw
elseif size(img,3)==3 && ...
        (all(elem(img(:,:,1)) == elem(img(:,:,2))) &&...
        all(elem(img(:,:,1)) == elem(img(:,:,3)))) && ...
        numel(unique(img(:)))==2
      % Assign image to bw image
    setappdata(handles.ArcasGui_figure, 'bw_img',logical(img(:,:,1)));
    % Convert image
    setappdata(handles.ArcasGui_figure, 'img',img);
    % Enable trial buttons
    trial_params_changed(handles);
else
    % Clear binary image
    setappdata(handles.ArcasGui_figure, 'bw_img',[]);
    % Store image
    setappdata(handles.ArcasGui_figure, 'img',img);
    disable_trial_buttons(handles);
end

% Enable thresholding checkboxes
set(handles.redThreshold_checkbox,'Enable','on');
set(handles.greenThreshold_checkbox,'Enable','on');
set(handles.blueThreshold_checkbox,'Enable','on');

guidata(hObject, handles);

function traitTot_edit_Callback(hObject, eventdata, handles)
% hObject    handle to traitTot_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of traitTot_edit as text
%        str2double(get(hObject,'String')) returns contents of traitTot_edit as a double
trial_params_changed(handles);

% --- Executes during object creation, after setting all properties.
function traitTot_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to traitTot_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cellPerTrial_edit_Callback(hObject, eventdata, handles)
% hObject    handle to cellPerTrial_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cellPerTrial_edit as text
%        str2double(get(hObject,'String')) returns contents of cellPerTrial_edit as a double
trial_params_changed(handles);

% --- Executes during object creation, after setting all properties.
function cellPerTrial_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cellPerTrial_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in runSimulation_pushbutton.
function runSimulation_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to runSimulation_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Set Cancel button visible to use can stop simulation
set(handles.cancel_togglebutton,'Visible', 'on');
set(handles.cancel_togglebutton,'Value', 0);
c = onCleanup(@(x) set(handles.cancel_togglebutton,'Visible', 'off'));
drawnow;

% Get total trials for simulation, events per trial, and diameter of agent
tot_trials = str2double(get(handles.traitTot_edit, 'String'));
events_per_trial = str2double(get(handles.cellPerTrial_edit,'String'));
event_diam = str2double(get(handles.cellPixelDiam_edit,'String'));
umppix = str2double(get(handles.umppix_edit,'String')); 
% Get background img
bw_img = getappdata(handles.ArcasGui_figure, 'bw_img');
assert(~isempty(bw_img), 'bw_img not found.')

% Get estimated time and max number of trials that can be run at once
max_trials_per_cycle = getappdata(handles.ArcasGui_figure,'max_trials_per_cycle');

[sim_hit_mean, sim_hit_std, event_vasc_dists, ~,trials_hit_rate] = ...
    ArcasGui_monteCarloSim_Driver(bw_img, event_diam, umppix, ...
    tot_trials, events_per_trial, max_trials_per_cycle, handles);
if isempty(sim_hit_mean) || isempty(sim_hit_std);  
    set(handles.cancel_togglebutton,'Value', 0);
    set(handles.cancel_togglebutton,'Visible', 'off')
    return; 
end
  
% figure;  
% histogram(trials_hit_rate*25,0:events_per_trial,'Normalization','probability')
% xlabel('Number of Cells Colocalizing')
% ylabel('Probability')
% beautifyAxis(gca)
% keyboard
% Save cell_vasc_dist
setappdata(handles.ArcasGui_figure,'trials_hite_rate',trials_hit_rate);
setappdata(handles.ArcasGui_figure,'cell_vasc_dists',event_vasc_dists);
set(handles.showHistogram_pushbutton, 'Enable','on');
set(handles.plotConvergence_pushbutton, 'Enable','on');

% Display results in simple textbox
set(handles.totTrials_text,'String',sprintf('%.f Trials', tot_trials));
set(handles.hitMean_text,'String',sprintf('Mean Colocalization Fraction           %10.4f', sim_hit_mean));
set(handles.hitSTD_text,'String', sprintf('STD Colocalization Fraction             %10.4f', sim_hit_std));



function cellPixelDiam_edit_Callback(hObject, eventdata, handles)
% hObject    handle to cellPixelDiam_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cellPixelDiam_edit as text
%        str2double(get(hObject,'String')) returns contents of cellPixelDiam_edit as a double
trial_params_changed(handles);

% --- Executes during object creation, after setting all properties.
function cellPixelDiam_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cellPixelDiam_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in singleTrial_pushbutton.
function singleTrial_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to singleTrial_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
events_per_trial = str2double(get(handles.cellPerTrial_edit,'String'));
event_diam = str2double(get(handles.cellPixelDiam_edit,'String'));
tot_trials = str2double(get(handles.traitTot_edit, 'String'));
umppix = str2double(get(handles.umppix_edit,'String')); 

% pixel diameter of agent
event_pix_diam = ceil(event_diam/umppix);

% Get input image
bw_img = getappdata(handles.ArcasGui_figure, 'bw_img'); 

% Run single trial, estimate memory and time
[comp_img, max_trials_per_cycle, minutes_for_simulation, bytes_per_trial] = ... 
    ArcasGui_RunSingleTrial(bw_img,events_per_trial,tot_trials, event_diam,umppix);

% Display estimated run time
set(handles.timeEstimate_text, 'String', sprintf('%s: %.2f Minutes',...
    'Estimate Simulation Run Time ',minutes_for_simulation));
% Display max concurrent trials
setappdata(handles.ArcasGui_figure, 'bytes_per_trial',bytes_per_trial);
set(handles.maxTrials_text, 'String', sprintf('%s: %.f ',' Max Concurrent Trials',...
    max_trials_per_cycle));
setappdata(handles.ArcasGui_figure,'max_trials_per_cycle',max_trials_per_cycle);

% Display Example trials
imshow(comp_img,'Parent', handles.img_axes);
% keyboard

% Enable simulation button
set(handles.runSimulation_pushbutton,'Enable','on');
set(handles.totTrials_text,'String','')
set(handles.hitMean_text,'String','')
set(handles.hitSTD_text,'String','')
% keyboard

function disable_trial_buttons(handles)
set(handles.singleTrial_pushbutton,'Enable','off');
set(handles.runSimulation_pushbutton,'Enable','off');
set(handles.showHistogram_pushbutton, 'Enable','off');
set(handles.plotConvergence_pushbutton, 'Enable','off');
set(handles.cancel_togglebutton,'Value', 0);
set(handles.cancel_togglebutton,'Visible', 'off');

set(handles.binaryNotification_text,'String','Image NOT Binary');
set(handles.binaryNotification_text,'ForegroundColor',[1 0 0]);
set(handles.threshold_uipanel,'BackgroundColor',[1 0 0]);


function trial_params_changed(handles)
% Enable/disable buttons if trial parameters changed (single trial must be
% rerun
set(handles.singleTrial_pushbutton,'Enable','on');
set(handles.runSimulation_pushbutton,'Enable','off');
set(handles.showHistogram_pushbutton, 'Enable','off');
set(handles.plotConvergence_pushbutton, 'Enable','off');

set(handles.cancel_togglebutton,'Value', 0);
set(handles.cancel_togglebutton,'Visible', 'off');

set(handles.binaryNotification_text,'String','Image Binary');
set(handles.binaryNotification_text,'ForegroundColor',[0 1 0]);
set(handles.threshold_uipanel,'BackgroundColor',[60 60 60]./255);


% --- Executes on button press in showHistogram_pushbutton.
function showHistogram_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to showHistogram_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cell_vasc_dists = getappdata(handles.ArcasGui_figure,'cell_vasc_dists');
if isempty(cell_vasc_dists); return; end
%plot data
figure; 
hist(cell_vasc_dists)
title('Histogram of Minimum Pixel Distance of Each Cell to Image Foreground')
xlabel('Minimum Distance'); ylabel('Frequency'); grid on;
y = ylim; x = xlim; axis([0 x(2) y(1) y(2)]);



function umppix_edit_Callback(hObject, eventdata, handles)
% hObject    handle to umppix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of umppix_edit as text
%        str2double(get(hObject,'String')) returns contents of umppix_edit as a double


% --- Executes during object creation, after setting all properties.
function umppix_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to umppix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cancel_pushbutton.
function cancel_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cancel_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cancel_togglebutton.
function cancel_togglebutton_Callback(hObject, eventdata, handles)
% hObject    handle to cancel_togglebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cancel_togglebutton


% --- Executes on slider movement.
function redThreshold_slider_Callback(hObject, eventdata, handles)
% hObject    handle to redThreshold_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function redThreshold_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to redThreshold_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in redThreshold_checkbox.
function redThreshold_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to redThreshold_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setThresholdButtons(hObject, [],handles);


% --- Executes on button press in greenThreshold_checkbox.
function greenThreshold_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to greenThreshold_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setThresholdButtons(hObject, [],handles);

% --- Executes on button press in blueThreshold_checkbox.
function blueThreshold_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to blueThreshold_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setThresholdButtons(hObject, [],handles);

% --- Executes on slider movement.
function greenThreshold_slider_Callback(hObject, eventdata, handles)
% hObject    handle to greenThreshold_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function greenThreshold_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to greenThreshold_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function blueThreshold_slider_Callback(hObject, eventdata, handles)
% hObject    handle to blueThreshold_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function blueThreshold_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blueThreshold_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in redThreshold_togglebutton.
function redThreshold_togglebutton_Callback(hObject, eventdata, handles)
% hObject    handle to redThreshold_togglebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject, 'value')==1
   set(hObject,'String', '<'); 
else
    
   set(hObject,'String', '>')
end




function Threshold_togglebutton_Callback(hObject, eventdata, handles)
if get(hObject, 'value')==1; set(hObject,'String', '<'); 
else set(hObject,'String', '>')
end
% initiate slide callback
threshold_sliders_Callback(hObject, eventdata, handles);






function setThresholdButtons(hObject, eventdata, handles)
name = get(hObject,'Tag');
channel_color = regexp(name, '(.*?)Thres','once','tokens');
if get(hObject,'Value')==1; val = 'on'; else val='inactive'; end
set(handles.([channel_color{1} 'Threshold_slider']),'Enable', val);
set(handles.([channel_color{1} 'Threshold_togglebutton']),'Enable', val);
if ~isempty(getappdata(handles.ArcasGui_figure,'img'))
    threshold_sliders_Callback(hObject, eventdata, handles);
end

function threshold_sliders_Callback(hObject, eventdata, handles)
img = getappdata(handles.ArcasGui_figure,'img');
[nr, nc, ~] = size(img);

use_red = get(handles.redThreshold_checkbox, 'Value');
use_green = get(handles.greenThreshold_checkbox, 'Value');
use_blue = get(handles.blueThreshold_checkbox, 'Value');

show_red = get(handles.redVeiw_togglebutton,'Value');
show_green = get(handles.greenVeiw_togglebutton,'Value');
show_blue = get(handles.blueVeiw_togglebutton,'Value');

if use_red && show_red
    red_bw = eval(['img(:,:,1) ' get(handles.redThreshold_togglebutton,'String') ...
        sprintf(' %.2f', get(handles.redThreshold_slider,'Value')*intmax(class(img)))]);
else red_bw = false([nr nc]); 
end

if use_green && show_green
    green_bw = eval(['img(:,:,2) ' get(handles.greenThreshold_togglebutton,'String') ...
        sprintf('%.2f', get(handles.greenThreshold_slider,'Value')*intmax(class(img)))]);
else green_bw = false([nr nc]); 
end

if use_blue && show_blue
    blue_bw = eval(['img(:,:,3) ' get(handles.blueThreshold_togglebutton,'String') ...
        sprintf('%.2f', get(handles.blueThreshold_slider,'Value')*intmax(class(img)))]);
else blue_bw = false([nr nc]); 
end

bw_img = bwareaopen(red_bw | green_bw | blue_bw, ...
    str2double(get(handles.minPixelArea_edit,'String')));
setappdata(handles.ArcasGui_figure,'bw_img', bw_img);


img(cat(3, bw_img,bw_img,bw_img))=intmax(class(img));

if any([use_red use_blue use_green])
    trial_params_changed(handles);
else
    disable_trial_buttons(handles);
end

% Display image with visible channels
channelVeiw_togglebutton_Callback(hObject, eventdata, handles);




function minPixelArea_edit_Callback(hObject, eventdata, handles)
% hObject    handle to minPixelArea_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minPixelArea_edit as text
%        str2double(get(hObject,'String')) returns contents of minPixelArea_edit as a double


% --- Executes during object creation, after setting all properties.
function minPixelArea_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minPixelArea_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function channelVeiw_togglebutton_Callback(hObject, eventdata, handles)

show_red = get(handles.redVeiw_togglebutton,'Value');
show_green = get(handles.greenVeiw_togglebutton,'Value');
show_blue = get(handles.blueVeiw_togglebutton,'Value');

% Get bw_img and img
bw_img = getappdata(handles.ArcasGui_figure,'bw_img');
img = getappdata(handles.ArcasGui_figure,'img');
img(:,:, ~logical([show_red show_green show_blue])) = 0;
img(cat(3, bw_img,bw_img,bw_img))=intmax(class(img));

imshow(img, 'Parent', handles.img_axes);


% --- Executes on button press in redVeiw_togglebutton.
function redVeiw_togglebutton_Callback(hObject, eventdata, handles)
% hObject    handle to redVeiw_togglebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of redVeiw_togglebutton


% --- Executes on button press in greenVeiw_togglebutton.
function greenVeiw_togglebutton_Callback(hObject, eventdata, handles)
% hObject    handle to greenVeiw_togglebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of greenVeiw_togglebutton


% --- Executes on button press in blueVeiw_togglebutton.
function blueVeiw_togglebutton_Callback(hObject, eventdata, handles)
% hObject    handle to blueVeiw_togglebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of blueVeiw_togglebutton


% --------------------------------------------------------------------
function export_menu_Callback(hObject, eventdata, handles)
% hObject    handle to export_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function exportBinary_submenu_Callback(hObject, eventdata, handles)
% hObject    handle to exportBinary_submenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
bw_img = getappdata(handles.ArcasGui_figure, 'bw_img');

img_path = getappdata(handles.ArcasGui_figure, 'img_path');

[pathstr,name,ext] = fileparts(img_path);
bw_name = [name '_bw.tif'];
% [bw_name, pathstr] = uiputfile('*.tif','Save Binary Image', bw_name);

imwrite(bw_img, [pathstr '/' bw_name]);


% --------------------------------------------------------------------
function exportCurrentImage_submenu_Callback(hObject, eventdata, handles)
% hObject    handle to exportCurrentImage_submenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img = getimage(handles.img_axes);

img_path = getappdata(handles.ArcasGui_figure, 'img_path');

[pathstr,name,ext] = fileparts(img_path);
img_name = [name '_sc.tif'];

[img_name, pathstr] = uiputfile('*.tif','Save Current Image', img_name);

imwrite(img, [pathstr '/' img_name]);


% --- Executes on button press in plotConvergence_pushbutton.
function plotConvergence_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to plotConvergence_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get list of hit rate for each trial
trials_hit_rate = getappdata(handles.ArcasGui_figure,'trials_hite_rate')';

% Calculate cumulative Mean
cum_mean = cumsum(trials_hit_rate)./(1:numel(trials_hit_rate));

% Calculate cumulative STD
x = trials_hit_rate;
n=numel(trials_hit_rate);
cum_std = sqrt((cumsum(x.^2) - (1:n).*(cumsum(x)./(1:n)).^2)./(0:(n-1)));
cum_std(1)=0;

figure
hAx = plotyy(1:numel(trials_hit_rate),cum_mean,1:numel(trials_hit_rate),cum_std);
% title('Convergence of Mean and STD across Trials')

xlabel('Trials Completed')
ylabel(hAx(1),'Cumulative Mean') % left y-axis
ylabel(hAx(2),'Cumulative STD') % right y-axis
grid on
beautifyAxis(gca)
set(gcf,'Position', [100 100 260 260])
keyboard
