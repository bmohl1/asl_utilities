function varargout = preproc_inputvars_asl_GUI(varargin)
%preproc_inputvars_asl_GUI M-file for preproc_inputvars_asl_GUI.fig
%      preproc_inputvars_asl_GUI, by itself, creates a new preproc_inputvars_asl_GUI or raises the existing
%      singleton*.
%
%      H = preproc_inputvars_asl_GUI returns the handle to a new preproc_inputvars_asl_GUI or the handle to
%      the existing singleton*.
%
%      preproc_inputvars_asl_GUI('Property','Value',...) creates a new preproc_inputvars_asl_GUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to preproc_inputvars_asl_GUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      preproc_inputvars_asl_GUI('CALLBACK') and preproc_inputvars_asl_GUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in preproc_inputvars_asl_GUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help preproc_inputvars_asl_GUI

% Last Modified by GUIDE v2.5 23-Oct-2017 08:39:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @preproc_inputvars_asl_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @preproc_inputvars_asl_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before preproc_inputvars_asl_GUI is made visible.
function preproc_inputvars_asl_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for preproc_inputvars_asl_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes preproc_inputvars_asl_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = preproc_inputvars_asl_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
specialTempflag = get(handles.checkbox_specialTemp,'Value');
saveQCflag = get(handles.checkbox_saveQC,'Value');
gmflag = get(handles.check_gmMasking,'Value');
redosegflag = get(handles.checkbox_redo_segment,'Value');
preprocflag = get(handles.checkbox_ignore,'Value');
cancelflag = get(handles.pushbutton_cancel,'userdata');

varargout{1} = specialTempflag;
varargout{2} = saveQCflag;
varargout{3} = gmflag; 
varargout{5} = redosegflag;
varargout{4} = preprocflag;
varargout{6} = cancelflag;

delete(handles.figure1);


% --- Executes on button press in push_ok.
function push_ok_Callback(hObject, eventdata, handles)
% hObject    handle to push_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(handles.figure1);
else
    % The GUI is no longer waiting, just close it
    delete(handles.figure1);
end

% --- Executes during object creation, after setting all properties.
function fc_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fc_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in checkbox_specialTemp.
function checkbox_specialTemp_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_specialTemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_specialTemp

handles.specialTemp = get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in checkbox_saveQC.
function checkbox_saveQC_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_saveQC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_saveQC
handles.saveQC = get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in check_gm.
function check_gmMasking_Callback(hObject, eventdata, handles)
% hObject    handle to check_gm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_gm
handles.gm = get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in check_discard.
%function check_unwarp_Callback(hObject, eventdata, handles)
% hObject    handle to check_discard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_discard
%handles.unwarp = get(hObject,'Value');
%guidata(hObject, handles);


% --- Executes on button press in checkbox_ignore.
function checkbox_ignore_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ignore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ignore
handles.ignorePreproc = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 set(handles.pushbutton_cancel,'userdata',1); % only changes to 1, if the function is invoked with button press
if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, use UIRESUME
    uiresume(handles.figure1);
else
    % The GUI is no longer waiting, just close it
    delete(handles.figure1);
end

% --- Executes on button press in checkbox_redo_segment.
function checkbox_redo_segment_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_redo_segment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_redo_segment
handles.redoseg = get(hObject,'Value');
guidata(hObject, handles);
