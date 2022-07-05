    function varargout = beautify(varargin)
% BEAUTIFY MATLAB code for beautify.fig
%      BEAUTIFY, by itself, creates a new BEAUTIFY or raises the existing
%      singleton*.
%
%      H = BEAUTIFY returns the handle to a new BEAUTIFY or the handle to
%      the existing singleton*.
%
%      BEAUTIFY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BEAUTIFY.M with the given input arguments.
%
%      BEAUTIFY('Property','Value',...) creates a new BEAUTIFY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before beautify_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to beautify_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help beautify

% Last Modified by GUIDE v2.5 23-Jun-2022 22:29:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @beautify_OpeningFcn, ...
    'gui_OutputFcn',  @beautify_OutputFcn, ...
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


% --- Executes just before beautify is made visible.
function beautify_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to beautify (see VARARGIN)

% Choose default command line output for beautify
handles.output = hObject;
set(handles.xpt,'String',1);
set(handles.xslide,'Min',0);
set(handles.ypt,'String',1);
set(handles.yslide,'Min',0);
set(handles.zpt,'String',1);
set(handles.zslide,'Min',0);
set(handles.tpt,'String',1);
set(handles.tslide,'Min',0);
set(handles.tfpt,'String',1);
set(handles.tfslide,'Min',0);
set(handles.use_chop,'Value',1);
% load_Callback(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes beautify wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = beautify_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function xmin_Callback(hObject, eventdata, handles)
% hObject    handle to xmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xmin as text
%        str2double(get(hObject,'String')) returns contents of xmin as a double


% --- Executes during object creation, after setting all properties.
function xmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xmax_Callback(hObject, eventdata, handles)
% hObject    handle to xmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xmax as text
%        str2double(get(hObject,'String')) returns contents of xmax as a double


% --- Executes during object creation, after setting all properties.
function xmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in max_range.
function max_range_Callback(hObject, eventdata, handles)
% hObject    handle to max_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.plotmod.Value == 1
    X = evalin('base','Xnew');
else
    X = evalin('base','X');
end
S = size(X);
%Sets min ranges to 1
set(handles.xmin,'String',1);
set(handles.ymin,'String',1);
set(handles.zmin,'String',1);
set(handles.tmin,'String',1);
set(handles.cmin,'String',min(X(:)));
%set(handles.cmin,'String','-1');
%Sets max ranges
set(handles.xmax,'String',S(1));
set(handles.ymax,'String',S(2));
set(handles.zmax,'String',S(3));
if length(S) >= 4
    set(handles.tmax,'String',S(4));
else
    set(handles.tmax,'String',1);
end
set(handles.cmax,'String',max(X(:)));
%set(handles.cmax,'String',1);



% Hint: get(hObject,'Value') returns toggle state of max_range


% --- Executes on button press in use_chop.
function use_chop_Callback(hObject, eventdata, handles)
% hObject    handle to use_chop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_chop


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
if handles.use_chop.Value
    X = evalin('base','X_c');
else
    X = evalin('base','Jrecon');
end
X = real(X);
assignin('base','X',X);
s = size(X);
% Writes number of points in each dimension
set(handles.xtext,'String',s(1));
set(handles.xslide,'Max',s(1))
set(handles.xslide,'SliderStep',[1/s(1) 5/s(1)]);
set(handles.ytext,'String',s(2));
set(handles.yslide,'Max',s(2))
set(handles.yslide,'SliderStep',[1/s(2) 5/s(2)]);
set(handles.ztext,'String',s(3));
set(handles.zslide,'Max',s(3))
set(handles.zslide,'SliderStep',[1/s(3) 5/s(3)]);
if length(s) == 4
    set(handles.ttext,'String',s(4));
    set(handles.tslide,'Max',s(4))
    set(handles.tslide,'SliderStep',[1/s(4) 5/s(4)]);
else
    set(handles.ttext,'String',1);
end


% Sets ranges for each dimension
max_range_Callback(hObject, eventdata, handles)


% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in plotorig.
function plotorig_Callback(hObject, eventdata, handles)
% hObject    handle to plotorig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.plotmod.Value
    X = evalin('base','Xnew');
else
    % H = handles.varselect.String{handles.varselect.Value};
    vname = handles.var_name.String;
    X = evalin('base',vname);
    M = handles.timeselect.String{handles.timeselect.Value};
%     switch H
%        case 'Pressure'
%            X = evalin('base','pressure');
           xpt = str2double(handles.xpt.String);
           ypt = str2double(handles.ypt.String);
           tfpt = str2double(handles.tfpt.String);
           tpt = str2double(handles.tpt.String);
           zpt = str2double(handles.zpt.String);
%        case 'Current'
%            X = evalin('base','I');
%        case 'AE'
%            X = evalin('base','AE');
%        case 'Lead'
%            X = evalin('base','L');
%        case 'Medium'
%            X = evalin('base','medium');
%        case 'Intensity'
%            X = evalin('base','sensor_data.p2');
%    end
end
N = get(handles.plot2d,'Value');
c = [str2double(get(handles.cmin,'String')) str2double(get(handles.cmax,'String'))];

if N == 1 %XZ
    p1 = str2double(get(handles.ypt,'String'));
    x = str2double(get(handles.xmin,'String')):str2double(get(handles.xmax,'String'));
    y = str2double(get(handles.zmin,'String')):str2double(get(handles.zmax,'String'));
    if handles.p_max.Value
        X2 = squeeze(X(x,p1,y,:));
        for i = 1:size(X2,1)
            for j = 1:size(X2,2)
                imag(i,j) = max(abs(X2(i,j,:)));
            end
        end
    else
        p2 = str2double(get(handles.tfpt,'String'));
        imag = squeeze(X(x,p1,y,p2))';
    end
elseif N == 2 %YZ
    p1 = str2double(get(handles.xpt,'String'));
    p2 = str2double(get(handles.tfpt,'String'));
    x = str2double(get(handles.ymin,'String')):str2double(get(handles.ymax,'String'));
    y = str2double(get(handles.zmin,'String')):str2double(get(handles.zmax,'String'));
    imag = squeeze(X(p1,x,y,p2))';
elseif N == 3 %XY
    p1 = zpt; 
    p2 = tfpt;
    p3 = tpt;
    x = str2double(get(handles.xmin,'String')):str2double(get(handles.xmax,'String'));
    y = str2double(get(handles.ymin,'String')):str2double(get(handles.ymax,'String'));
    if handles.p_max.Value
        if ndims(X) == 3
            X2 = squeeze(X(x,y,:));
        elseif ndims(X) == 4
            X2 = squeeze(X(x,y,p1,:));
        end
        for i = 1:size(X2,1)
            for j = 1:size(X2,2)
                imag(i,j) = max(abs(X2(i,j,:)));
            end
            multiWaitbar('Solving Max Signal',i/size(X2,1));
        end
        multiWaitbar('CLOSEALL')
    else
        if ndims(X) == 3    
            p4 = max([p1 p2 p3]);
            imag = squeeze(X(x,y,p4));
        elseif ndims(X) == 4
            imag = squeeze(X(x,y,p1,p2));
        else

            imag = squeeze(X(x,y,p1,p2,p3));

        end
    end
    %%%%%%%%%%%%% INDIVIDUAL IMAGE %%%%%%%%%%%%
    %cmax = max(imag(:))*c(2);
    %cmin = min(imag(:))*c(1)*-1;

elseif N == 4 %TX
    p1 = str2double(get(handles.ypt,'String'));
    p2 = str2double(get(handles.zpt,'String'));
    x = str2double(get(handles.xmin,'String')):str2double(get(handles.xmax,'String'));
    y = str2double(get(handles.tfmin,'String')):str2double(get(handles.tfmax,'String'));
    imag = squeeze(X(x,p1,p2,y));
elseif N == 5 %TY
    p1 = str2double(get(handles.xpt,'String'));
    p2 = str2double(get(handles.zpt,'String'));
    x = str2double(get(handles.ymin,'String')):str2double(get(handles.ymax,'String'));
    y = str2double(get(handles.tfmin,'String')):str2double(get(handles.tfmax,'String'));
    imag = squeeze(X(p1,x,p2,y));
elseif N == 6 %TZ
    p1 = str2double(get(handles.xpt,'String'));
    p2 = str2double(get(handles.ypt,'String'));
    x = str2double(get(handles.zmin,'String')):str2double(get(handles.zmax,'String'));
    y = str2double(get(handles.tfmin,'String')):str2double(get(handles.tfmax,'String'));
    imag = squeeze(X(p1,p2,x,y));
end
hmap = get(handles.cmapmenu,'String');
h = hmap{get(handles.cmapmenu,'Value')};
if get(handles.cmapmenu,'Value') == 4
    clear h
    h = hotcoldDB;
elseif get(handles.cmapmenu,'Value') == 3
    clear h
    h = hotcold;
end
if handles.ext_fig.Value
    figure(21);
else
    axes(handles.axes1)
end
if handles.op_transpose.Value
    imag = imag';
end
if handles.op_env.Value
    imag = envelope(imag')';
end
if handles.bband.Value
%     template = str2double(handles.bb_temp.String);
%     T = imag(:,template);
%     L = size(imag,1);
%     for i = 1:size(imag,2)
%         imag2(:,i) = xcorr(imag(:,i),T);
%     end
%     clear imag
%     for i = 1:size(imag2,2)
%         imag(:,i) = interp1(linspace(0,1,size(imag2,1)),imag2(:,i),linspace(0,1,L));
%     end
fc = str2double(handles.bb_temp.String);
dy = str2double(handles.dy.String);
fs = str2double(handles.fs.String);
L1 = size(imag,1);
T = linspace(0,L1/fs,L1);
imag = ae_demod3(imag,T,fc); 


end
if handles.filter.Value
    m = [str2double(handles.med_x.String) str2double(handles.med_y.String) str2double(handles.med_z.String)];
    n = [str2double(handles.mean_x.String) str2double(handles.mean_y.String) str2double(handles.mean_z.String)];
    o = [str2double(handles.int_x.String) str2double(handles.int_y.String) str2double(handles.int_z.String)];
    switch N
        case 1
            m1 = m(1);
            m2 = m(3);
            n1 = n(1);
            n2 = n(3);
            o1 = o(1);
            o2 = o(3);
        case 2
            m1 = m(2);
            m2 = m(3);
            n1 = n(2);
            n2 = n(3);
            o1 = o(2);
            o2 = o(3);
        case 3
            m1 = m(1);
            m2 = m(2);
            n1 = n(1);
            n2 = n(2);
            o1 = o(1);
            o2 = o(2);
    end
    if m1+m2 > 2
    end
    if n1+n2 > 2
        imag = imgaussfilt(imag,(n1+n2)/2);
    end
end
if handles.apply_mods.Value
    m = [0 str2double(handles.med_x.String) str2double(handles.med_y.String) str2double(handles.med_z.String)];
    n = [0 str2double(handles.mean_x.String) str2double(handles.mean_y.String) str2double(handles.mean_z.String)];
    o = [0 str2double(handles.int_x.String) str2double(handles.int_y.String) str2double(handles.int_z.String)];
    if sum(m) > 3
        m(1) = 1;
    end
    if sum(n) > 3
        n(1) = 1;
    end
    if sum(o) > 3
        o(1) = 1;
    end
    param.window = 'beautify';
    p = [1 1 1];
    imag = filts3D(imag,n,o,m,p,param);
end
if handles.op_invert.Value
    imag = imag*-1;
end
if handles.op_transpose.Value
    imag = imag';
end

    %%%%%%%%%% GLOBAL MIN MAX %%%%%%%%%%%%%%
    if handles.dbc.Value
        cmax = c(2);
        cmin = c(1);
        if strcmp(handles.cmapmenu.String{handles.cmapmenu.Value},'hotcoldDB')
            S = sign(imag);
            imag = real(20*log10(imag./max(imag(:))));
            imag = S.*imag;
        else
            imag = real(20*log10(imag./max(imag(:))));
        end
    else
        cmax = max(imag(:))*c(2);
        cmin = min(imag(:))*c(1)*-1;
    end
    c = [cmin cmax];
if handles.convert.Value
    xcon = str2double(handles.dx.String);
    ycon = str2double(handles.dy.String);
    xax = y.*xcon;
    yax = x.*ycon;
else
    xax = y;
    yax = x;
end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isnan(c)
    if c(1) == c(2)
        imagesc(xax,yax,imag,'ButtonDownFcn',{@Plot_FWHM,handles});
    else
        imagesc(xax,yax,imag,'ButtonDownFcn',{@Plot_FWHM,handles},c);
    end
else
    imagesc(xax,yax,imag,'ButtonDownFcn',{@Plot_FWHM,handles});
end
colormap(h)
% c = [str2double(get(handles.cmin,'String')) str2double(get(handles.cmax,'String'))];
% if isnan(c)
% else
%     caxis(c);
% end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function xpt_Callback(hObject, eventdata, handles)
% hObject    handle to xpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xpt as text
%        str2double(get(hObject,'String')) returns contents of xpt as a double


% --- Executes during object creation, after setting all properties.
function xpt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ymin_Callback(hObject, eventdata, handles)
% hObject    handle to ymin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ymin as text
%        str2double(get(hObject,'String')) returns contents of ymin as a double


% --- Executes during object creation, after setting all properties.
function ymin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ymin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ymax_Callback(hObject, eventdata, handles)
% hObject    handle to ymax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ymax as text
%        str2double(get(hObject,'String')) returns contents of ymax as a double


% --- Executes during object creation, after setting all properties.
function ymax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ymax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ypt_Callback(hObject, eventdata, handles)
% hObject    handle to ypt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ypt as text
%        str2double(get(hObject,'String')) returns contents of ypt as a double


% --- Executes during object creation, after setting all properties.
function ypt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ypt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zmin_Callback(hObject, eventdata, handles)
% hObject    handle to zmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zmin as text
%        str2double(get(hObject,'String')) returns contents of zmin as a double


% --- Executes during object creation, after setting all properties.
function zmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zmax_Callback(hObject, eventdata, handles)
% hObject    handle to zmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zmax as text
%        str2double(get(hObject,'String')) returns contents of zmax as a double


% --- Executes during object creation, after setting all properties.
function zmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zpt_Callback(hObject, eventdata, handles)
% hObject    handle to zpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zpt as text
%        str2double(get(hObject,'String')) returns contents of zpt as a double


% --- Executes during object creation, after setting all properties.
function zpt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tmin_Callback(hObject, eventdata, handles)
% hObject    handle to tmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tmin as text
%        str2double(get(hObject,'String')) returns contents of tmin as a double


% --- Executes during object creation, after setting all properties.
function tmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tmax_Callback(hObject, eventdata, handles)
% hObject    handle to tmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tmax as text
%        str2double(get(hObject,'String')) returns contents of tmax as a double


% --- Executes during object creation, after setting all properties.
function tmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tpt_Callback(hObject, eventdata, handles)
% hObject    handle to tpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tpt as text
%        str2double(get(hObject,'String')) returns contents of tpt as a double


% --- Executes during object creation, after setting all properties.
function tpt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in plot2d.
function plot2d_Callback(hObject, eventdata, handles)
% hObject    handle to plot2d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plot2d contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plot2d


% --- Executes during object creation, after setting all properties.
function plot2d_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot2d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in plot3d.
function plot3d_Callback(hObject, eventdata, handles)
% hObject    handle to plot3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plot3d contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plot3d


% --- Executes during object creation, after setting all properties.
function plot3d_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cmin_Callback(hObject, eventdata, handles)
% hObject    handle to cmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cmin as text
%        str2double(get(hObject,'String')) returns contents of cmin as a double


% --- Executes during object creation, after setting all properties.
function cmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cmax_Callback(hObject, eventdata, handles)
% hObject    handle to cmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cmax as text
%        str2double(get(hObject,'String')) returns contents of cmax as a double


% --- Executes during object creation, after setting all properties.
function cmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in cmapmenu.
function cmapmenu_Callback(hObject, eventdata, handles)
% hObject    handle to cmapmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cmapmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cmapmenu


% --- Executes during object creation, after setting all properties.
function cmapmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cmapmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in moviebutton.
function moviebutton_Callback(hObject, eventdata, handles)
% hObject    handle to moviebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.plotmod.Value
    X = evalin('base','Xnew');
else H = handles.varselect.String{handles.varselect.Value};
%     switch H
%         case 'Pressure'
%             X = evalin('base','pressure');
%         case 'Current'
%             X = evalin('base','I');
%         case 'AE'
%             X = evalin('base','AE');
%         case 'Lead'
%             X = evalin('base','L');
%         case 'Medium'
%             X = evalin('base','medium');
%         case 'Intensity'
%             X = evalin('base','sensor_data.p2');
%     end
    vname = handles.var_name.String;
    X = evalin('base',vname);
end
N = get(handles.plot3d,'Value');

% if ndims(X) == 3
%     X = permute(X,[1,2,4,3]);
% end
if N == 1 %XZt
    p = str2double(get(handles.ypt,'String'));
    x = str2double(get(handles.xmin,'String')):str2double(get(handles.xmax,'String'));
    y = str2double(get(handles.zmin,'String')):str2double(get(handles.zmax,'String'));
    t = str2double(get(handles.tfmin,'String')):str2double(get(handles.tfmax,'String'));
     tstart = str2double(get(handles.tfmin,'String'));
    imag = squeeze(X(x,p,y,:));
elseif N == 2 %YZt
    p = str2double(get(handles.xpt,'String'));
    x = str2double(get(handles.ymin,'String')):str2double(get(handles.ymax,'String'));
    y = str2double(get(handles.zmin,'String')):str2double(get(handles.zmax,'String'));
    t = str2double(get(handles.tfmin,'String')):str2double(get(handles.tfmax,'String'));
     tstart = str2double(get(handles.tfmin,'String'));
    imag = squeeze(X(p,x,y,:));
elseif N == 3 %XYt
    p = str2double(get(handles.zpt,'String'));
    x = str2double(get(handles.xmin,'String')):str2double(get(handles.xmax,'String'));
    y = str2double(get(handles.ymin,'String')):str2double(get(handles.ymax,'String'));
    t = str2double(get(handles.tfmin,'String')):str2double(get(handles.tfmax,'String'));
     tstart = str2double(get(handles.tfmin,'String'));
    imag = squeeze(X(x,y,p,:));
elseif N == 4 %XYz
    p = str2double(get(handles.tfpt,'String'));
    x = str2double(get(handles.xmin,'String')):str2double(get(handles.xmax,'String'));
    y = str2double(get(handles.ymin,'String')):str2double(get(handles.ymax,'String'));
    t = str2double(get(handles.zmin,'String')):str2double(get(handles.zmax,'String'));
    tstart = str2double(get(handles.zmin,'String'));
    imag = squeeze(X(x,y,t,p));
elseif N == 5 %XZy
    p = str2double(get(handles.tfpt,'String'));
    x = str2double(get(handles.xmin,'String')):str2double(get(handles.xmax,'String'));
    y = str2double(get(handles.zmin,'String')):str2double(get(handles.zmax,'String'));
    t = str2double(get(handles.ymin,'String')):str2double(get(handles.ymax,'String'));
     tstart = str2double(get(handles.ymin,'String'));
    imag = squeeze(X(x,:,y,p));
    imag = permute(imag,[1 3 2]);
elseif N == 6 %YZx
    p = str2double(get(handles.tfpt,'String'));
    x = str2double(get(handles.ymin,'String')):str2double(get(handles.ymax,'String'));
    y = str2double(get(handles.zmin,'String')):str2double(get(handles.zmax,'String'));
    t = str2double(get(handles.xmin,'String')):str2double(get(handles.xmax,'String'));
     tstart = str2double(get(handles.xmin,'String'));
    imag = squeeze(X(:,x,y,p));
    imag = permute(imag,[2 3 1]);
end
c = [str2double(get(handles.cmin,'String')) str2double(get(handles.cmax,'String'))];
c2 = [c(1)*min(imag(:)) c(2)*max(imag(:))];
hmap = get(handles.cmapmenu,'String');
h = hmap{get(handles.cmapmenu,'Value')};
if get(handles.cmapmenu,'Value') == 4
    clear h
    h = hotcoldDB;
elseif get(handles.cmapmenu,'Value') == 3
    clear h
    h = hotcold;
end
if handles.savebox.Value
    vidwrite2(imag(:,:,:),handles,c,h);
else
    if handles.plotmod.Value
        axes(handles.axes3)
    else
        axes(handles.axes1)
    end
    if handles.op_env.Value
        for i = 1:size(imag,3)
            a(:,:,i) = envelope(imag(:,:,i));
            multiWaitbar('Enveloping',i/size(imag,3));
        end
        multiWaitbar('CLOSEALL');
        imag = a;
    end
    if handles.op_transpose.Value
        imag = permute(imag,[2,1,3]);
    end
    if handles.op_invert.Value
        imag = imag*-1;
    end
    for i = t-tstart+1
        if i == t(1)-tstart+1
            cmax = max(imag(:))*c(2);
            cmin = min(imag(:))*c(1)*-1;
            c = [cmin cmax];
            if ~isnan(c)
                imagesc(imag(:,:,i),c);
            else
            imagesc(imag(:,:,i));
            end
        else
            if handles.plotmod.Value
                handles.axes3.Children.CData  = imag(:,:,i);
            else
                handles.axes1.Children.CData  = imag(:,:,i);
            end
        end
        colormap(h);
        
        if isnan(c)
        else
            caxis(c);
        end
        title(i);
        %drawnow;
        pause(str2double(handles.pausems.String)/1e3);
    end
end



% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in holdmods.
function holdmods_Callback(hObject, eventdata, handles)
% hObject    handle to holdmods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of holdmods



function med_x_Callback(hObject, eventdata, handles)
% hObject    handle to med_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of med_x as text
%        str2double(get(hObject,'String')) returns contents of med_x as a double


% --- Executes during object creation, after setting all properties.
function med_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to med_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function med_y_Callback(hObject, eventdata, handles)
% hObject    handle to med_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of med_y as text
%        str2double(get(hObject,'String')) returns contents of med_y as a double


% --- Executes during object creation, after setting all properties.
function med_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to med_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function med_z_Callback(hObject, eventdata, handles)
% hObject    handle to med_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of med_z as text
%        str2double(get(hObject,'String')) returns contents of med_z as a double


% --- Executes during object creation, after setting all properties.
function med_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to med_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tmean_Callback(hObject, eventdata, handles)
% hObject    handle to tmean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tmean as text
%        str2double(get(hObject,'String')) returns contents of tmean as a double


% --- Executes during object creation, after setting all properties.
function tmean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tmean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ext_fig.
function ext_fig_Callback(hObject, eventdata, handles)
% hObject    handle to ext_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ext_fig


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.holdmods.Value
    X = evalin('base','Xnew');
else
    H = handles.varselect.String{handles.varselect.Value};
    switch H
        case 'Pressure'
            X = evalin('base','pressure');
        case 'Current'
            X = evalin('base','I');
        case 'AE'
            X = evalin('base','AE');
        case 'Lead'
            X = evalin('base','L');
        case 'Medium'
            X = evalin('base','medium');
        case 'Intensity'
            X = evalin('base','sensor_data.p2');
    end
end
x = [str2double(get(handles.xmin,'String')) str2double(get(handles.xmax,'String')) str2double(get(handles.xpt,'String'))];
y = [str2double(get(handles.ymin,'String')) str2double(get(handles.ymax,'String')) str2double(get(handles.ypt,'String'))];
z = [str2double(get(handles.zmin,'String')) str2double(get(handles.zmax,'String')) str2double(get(handles.zpt,'String'))];
t = [str2double(get(handles.tmin,'String')) str2double(get(handles.tmax,'String')) str2double(get(handles.tpt,'String'))];

%Xmod = X(x(1):x(2),y(1):y(2),z(1):z(2),t(1):t(2));

box = [str2double(get(handles.med_x,'String')) str2double(get(handles.med_y,'String')) str2double(get(handles.med_z,'String'))];

if handles.allt.Value == 1
    n = t(1):t(2);
else
    n = t(3);
end
for i = n
    Xtemp(:,:,:,i) = imboxfilt3(X(x(1):x(2),y(1):y(2),z(1):z(2),i),box);%imboxfilt3(Xmod(:,:,:,n),box);
end

X(x(1):x(2),y(1):y(2),z(1):z(2),n) = Xtemp(:,:,:,n);

assignin('base','Xnew',X);

if handles.holdmods.Value
    X = evalin('base','Xnew');
else
    X = evalin('base',handles.var_name.String);
end
x = [str2double(get(handles.xmin,'String')) str2double(get(handles.xmax,'String')) str2double(get(handles.xpt,'String'))];
y = [str2double(get(handles.ymin,'String')) str2double(get(handles.ymax,'String')) str2double(get(handles.ypt,'String'))];
z = [str2double(get(handles.zmin,'String')) str2double(get(handles.zmax,'String')) str2double(get(handles.zpt,'String'))];
t = [str2double(get(handles.tmin,'String')) str2double(get(handles.tmax,'String')) str2double(get(handles.tpt,'String'))];

m1 = [str2double(handles.mean_x.String) str2double(handles.mean_y.String) str2double(handles.mean_z.String)];
n1 = [str2double(handles.int_x.String) str2double(handles.int_y.String) str2double(handles.int_z.String) 0];% get(handles.squarify_box,'Value')];
o1 = [str2double(handles.med_x.String) str2double(handles.med_y.String) str2double(handles.med_z.String)];
if sum(m1) ~= 3
    m = [1 m1];
else
    m = [0 m1];
end
if sum(n1) ~= 3
    n = [1 n1];
else
    n = [0 n1];
end
if sum(o1) ~= 3
    o = [1 o1];
else
    o = [0 o1];
end
param.window = '1D Guass';
p = [1 1 1];
Xtemp = X(x(1):x(2),y(1):y(2),z(1):z(2),t(1):t(2));
X2 = filts3D(Xtemp,m,n,o,p,param);
%    X(x(1):round(x(2)*o(1)),y(1):round(y(2)*o(2)),z(1):round(z(2)*o(3)),t(1):t(2)) = X2;
assignin('base','Xnew',X2)


% --- Executes on button press in allt.
function allt_Callback(hObject, eventdata, handles)
% hObject    handle to allt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of allt


% --- Executes on button press in plotmod.
function plotmod_Callback(hObject, eventdata, handles)
% hObject    handle to plotmod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if hObject.Value
    set(handles.plot5d,'BackgroundColor',[0.7 0 0])
    set(handles.movie5d,'BackgroundColor',[0.7 0 0])
    set(handles.plotorig,'BackgroundColor',[0.7 0 0])
    set(handles.moviebutton,'BackgroundColor',[0.7 0 0])
    set(handles.plot5d,'ForegroundColor',[0.9 0.9 0.9])
    set(handles.movie5d,'ForegroundColor',[0.9 0.9 0.9])
    set(handles.plotorig,'ForegroundColor',[0.9 0.9 0.9])
    set(handles.moviebutton,'ForegroundColor',[0.9 0.9 0.9])
else
    set(handles.plot5d,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.movie5d,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.plotorig,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.moviebutton,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.plot5d,'ForegroundColor',[0 0 0])
    set(handles.movie5d,'ForegroundColor',[0 0 0])
    set(handles.plotorig,'ForegroundColor',[0 0 0])
    set(handles.moviebutton,'ForegroundColor',[0 0 0])
end
% Hint: get(hObject,'Value') returns toggle state of plotmod


% --- Executes on button press in clearc.
function clearc_Callback(hObject, eventdata, handles)
% hObject    handle to clearc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.cmax,'String',[]);
set(handles.cmin,'String',[]);
% Hint: get(hObject,'Value') returns toggle state of clearc


% --- Executes on button press in loadparams.
function loadparams_Callback(hObject, eventdata, handles)
% hObject    handle to loadparams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s = evalin('base','selection');
set(handles.xmin,'String',s.x(1));
set(handles.xmax,'String',s.x(2));
set(handles.xpt,'String',s.x(3));
set(handles.ymin,'String',s.y(1));
set(handles.ymax,'String',s.y(2));
set(handles.ypt,'String',s.y(3));
set(handles.zmin,'String',s.z(1));
set(handles.zmax,'String',s.z(2));
set(handles.zpt,'String',s.z(3));
set(handles.tmin,'String',s.t(1));
set(handles.tmax,'String',s.t(2));
set(handles.tpt,'String',s.t(3));
if isnan(s.c)
    s.c = [];
    set(handles.cmin,'String',s.c);
    set(handles.cmax,'String',s.c);
else
    set(handles.cmin,'String',s.c(1));
    set(handles.cmax,'String',s.c(2));
end



% --- Executes on button press in saveparams.
function saveparams_Callback(hObject, eventdata, handles)
% hObject    handle to saveparams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection.x = [str2double(get(handles.xmin,'String')) str2double(get(handles.xmax,'String')) str2double(get(handles.xpt,'String'))];
selection.y = [str2double(get(handles.ymin,'String')) str2double(get(handles.ymax,'String')) str2double(get(handles.ypt,'String'))];
selection.z = [str2double(get(handles.zmin,'String')) str2double(get(handles.zmax,'String')) str2double(get(handles.zpt,'String'))];
selection.t = [str2double(get(handles.tmin,'String')) str2double(get(handles.tmax,'String')) str2double(get(handles.tpt,'String'))];
selection.c = [str2double(get(handles.cmin,'String')) str2double(get(handles.cmax,'String'))];
assignin('base','selection',selection);



function pausems_Callback(hObject, eventdata, handles)
% hObject    handle to pausems (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pausems as text
%        str2double(get(hObject,'String')) returns contents of pausems as a double


% --- Executes during object creation, after setting all properties.
function pausems_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pausems (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.holdmods.Value
    X = evalin('base','Xnew');
else
 H = handles.varselect.String{handles.varselect.Value};
    switch H
        case 'Pressure'
            X = evalin('base','pressure');
        case 'Current'
            X = evalin('base','I');
        case 'AE'
            X = evalin('base','AE');
        case 'Lead'
            X = evalin('base','L');
        case 'Medium'
            X = evalin('base','medium');
        case 'Intensity'
            X = evalin('base','sensor_data.p2');
    end
end
x = [str2double(get(handles.xmin,'String')) str2double(get(handles.xmax,'String')) str2double(get(handles.xpt,'String'))];
y = [str2double(get(handles.ymin,'String')) str2double(get(handles.ymax,'String')) str2double(get(handles.ypt,'String'))];
z = [str2double(get(handles.zmin,'String')) str2double(get(handles.zmax,'String')) str2double(get(handles.zpt,'String'))];
t = [str2double(get(handles.tfmin,'String')) str2double(get(handles.tfmax,'String')) str2double(get(handles.tfpt,'String'))];
if handles.allt.Value
    Xnew = abs(hilbert(X(x(1):x(2),y(1):y(2),z(1):z(2),:)));
else
    Xnew = zeros(x(2)-x(1)+1,y(2)-y(1)+1,z(2)-z(1)+1,t(2)-t(1)+1);
    for i = y(1):y(2)
        for j = z(1):z(2)
            for k = t(1):t(2)
                Xnew(:,i,j,k) = abs(hilbert(X(x(1):x(2),i,j,k)));
            end
        end
        multiWaitbar('Enveloping',i/(y(2)-y(1)+1));
    end
    %Xnew = abs(X);
    multiWaitbar('CLOSE ALL');
end
assignin('base','Xnew',Xnew);



function mean_x_Callback(hObject, eventdata, handles)
% hObject    handle to mean_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mean_x as text
%        str2double(get(hObject,'String')) returns contents of mean_x as a double


% --- Executes during object creation, after setting all properties.
function mean_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mean_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mean_y_Callback(hObject, eventdata, handles)
% hObject    handle to mean_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mean_y as text
%        str2double(get(hObject,'String')) returns contents of mean_y as a double


% --- Executes during object creation, after setting all properties.
function mean_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mean_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mean_z_Callback(hObject, eventdata, handles)
% hObject    handle to mean_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mean_z as text
%        str2double(get(hObject,'String')) returns contents of mean_z as a double


% --- Executes during object creation, after setting all properties.
function mean_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mean_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function int_x_Callback(hObject, eventdata, handles)
% hObject    handle to int_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of int_x as text
%        str2double(get(hObject,'String')) returns contents of int_x as a double


% --- Executes during object creation, after setting all properties.
function int_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to int_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function int_y_Callback(hObject, eventdata, handles)
% hObject    handle to int_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of int_y as text
%        str2double(get(hObject,'String')) returns contents of int_y as a double


% --- Executes during object creation, after setting all properties.
function int_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to int_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function int_z_Callback(hObject, eventdata, handles)
% hObject    handle to int_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of int_z as text
%        str2double(get(hObject,'String')) returns contents of int_z as a double


% --- Executes during object creation, after setting all properties.
function int_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to int_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.plotmod.Value
    X = evalin('base','Xnew');
else
    X = evalin('base','X');
end
x = str2num(handles.xmin.String):str2num(handles.xmax.String);
y = str2num(handles.ymin.String):str2num(handles.ymax.String);
z = str2num(handles.zmin.String):str2num(handles.zmax.String);
t = str2num(handles.tmin.String):str2num(handles.tmax.String);
X = X(x,y,z,t);
assignin('base','X_c',X);


% --- Executes on slider movement.
function xslide_Callback(hObject, eventdata, handles)
% hObject    handle to xslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = round(get(hObject,'Value'));
t = [str2double(handles.xmin.String) str2double(handles.xmax.String)];
if val > t(2)
    val = t(2)-1;
end
if val < t(1)
    val = t(1)+1;
end
handles.xslide.Max = t(2);
handles.xslide.Min = t(1);
set(handles.xslide,'SliderStep',[1/t(2) 5/t(2)]);
set(handles.xpt,'String',val);
set(handles.xslide,'Value',val);
if handles.silder5d.Value
    plot5d_Callback(hObject, eventdata,handles);
else
    plotorig_Callback(hObject, eventdata, handles)
end


%handles.xslide.SliderStep = [1 10];

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function xslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function yslide_Callback(hObject, eventdata, handles)
% hObject    handle to yslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = round(get(hObject,'Value'));
t = [str2double(handles.ymin.String) str2double(handles.ymax.String)];
if val > t(2)
    val = t(2)-1;
end
if val < t(1)
    val = t(1)+1;
end
handles.yslide.Max = t(2);
handles.yslide.Min = t(1);
set(handles.yslide,'SliderStep',[1/t(2) 5/t(2)]);
set(handles.ypt,'String',val);
set(handles.yslide,'Value',val);
if handles.silder5d.Value
    plot5d_Callback(hObject, eventdata,handles);
else
    plotorig_Callback(hObject, eventdata, handles)
end
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function yslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function zslide_Callback(hObject, eventdata, handles)
% hObject    handle to zslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = round(get(hObject,'Value'));
t = [str2double(handles.zmin.String) str2double(handles.zmax.String)];
if val > t(2)
    val = t(2)-1;
end
if val < t(1)
    val = t(1)+1;
end
handles.zslide.Max = t(2);
handles.zslide.Min = t(1);
set(handles.zslide,'SliderStep',[1/t(2) 5/t(2)]);
set(handles.zpt,'String',val);
set(handles.zslide,'Value',val);
if handles.silder5d.Value
    plot5d_Callback(hObject, eventdata,handles);
else
    plotorig_Callback(hObject, eventdata, handles)
end
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function zslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function tslide_Callback(hObject, eventdata, handles)
% hObject    handle to tslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = round(get(hObject,'Value'));
t = [str2double(handles.tmin.String) str2double(handles.tmax.String)];
if val > t(2)
    val = t(2)-1;
end
if val < t(1)
    val = t(1)+1;
end
handles.tslide.Max = t(2);
handles.tslide.Min = t(1);
set(handles.tslide,'SliderStep',[1/t(2) 5/t(2)]);
set(handles.tpt,'String',val);
set(handles.tslide,'Value',val);
if handles.silder5d.Value
    plot5d_Callback(hObject, eventdata,handles);
else
    plotorig_Callback(hObject, eventdata, handles)
end
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function tslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in TodB.
function TodB_Callback(hObject, eventdata, handles)
% hObject    handle to TodB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.holdmods.Value
    X = evalin('base','Xnew');
else
    X = evalin('base','X');
end
S = sign(X);
D = abs(X);
B = 20*log10(D./max(D(:)));
C  = B.*S;
assignin('base','Xnew',C);


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.plotmod.Value
    X = evalin('base','Xnew');
else
    X = evalin('base','X');
end
d = size(X);
if ndims(X) == 3
    if handles.vol.Value
        d = [d 1];
    else
        d = [d(1) 1 d(2) d(3)];
    end
elseif ndims(X) == 2
    d = [d(1) 1 d(2) 1];
end
x = str2num(handles.xmin.String):str2num(handles.xmax.String);
y = str2num(handles.ymin.String):str2num(handles.ymax.String);
z = str2num(handles.zmin.String):str2num(handles.zmax.String);
t = str2num(handles.tmin.String):str2num(handles.tmax.String);
a.x = x;
a.y = y;
a.depth = z;
a.stime = t;
ax_c.x = 1:d(1);
ax_c.y = 1:d(2);
ax_c.depth = 1:d(3);
ax_c.stime = 1:d(4);
assignin('base','ax_c',a);


% --- Executes on button press in vol.
function vol_Callback(hObject, eventdata, handles)
% hObject    handle to vol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vol


% --- Executes on button press in savebox.
function savebox_Callback(hObject, eventdata, handles)
% hObject    handle to savebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of savebox



function savefolder_Callback(hObject, eventdata, handles)
% hObject    handle to savefolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of savefolder as text
%        str2double(get(hObject,'String')) returns contents of savefolder as a double


% --- Executes during object creation, after setting all properties.
function savefolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savefolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function savefigname_Callback(hObject, eventdata, handles)
% hObject    handle to savefigname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of savefigname as text
%        str2double(get(hObject,'String')) returns contents of savefigname as a double


% --- Executes during object creation, after setting all properties.
function savefigname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savefigname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in thisfold.
function thisfold_Callback(hObject, eventdata, handles)
% hObject    handle to thisfold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.savefolder,'String',pwd);



function framerate_Callback(hObject, eventdata, handles)
% hObject    handle to framerate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of framerate as text
%        str2double(get(hObject,'String')) returns contents of framerate as a double


% --- Executes during object creation, after setting all properties.
function framerate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to framerate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
X = evalin('base','X');
D = size(X);
if length(D) == 2
    D(3:4) = 1;
elseif length(D) == 3
    D(4) = 1;
end
m1 = [str2double(handles.mean_x.String) str2double(handles.mean_y.String) str2double(handles.mean_z.String)];
n1 = [str2double(handles.int_x.String) str2double(handles.int_y.String) str2double(handles.int_z.String) 0];% get(handles.squarify_box,'Value')];
o1 = [str2double(handles.med_x.String) str2double(handles.med_y.String) str2double(handles.med_z.String)];
if sum(m1) ~= 3
    m = [1 m1];
else
    m = [0 m1];
end
if sum(n1) ~= 3
    n = [1 n1];
else
    n = [0 n1];
end
if sum(o1) ~= 3
    o = [1 o1];
else
    o = [0 o1];
end
n2 = [0 n1];
X2 = filts3D(X,m,n2,o);
S = sign(X2);
X2 = X2.^2;
X2 = X2.*S;
X3 = filts3D(X2,m,n,o);

assignin('base','Xnew',X3);
max_range_Callback(hObject, eventdata, handles)
pushbutton11_Callback(hObject, eventdata, handles)
X = evalin('base','Xnew');
X4 = filts3D(X,m,n2,o);
assignin('base','Xnew',X4);
max_range_Callback(hObject, eventdata, handles)
a = 3;

% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    


% --- Executes on button press in shift.
function shift_Callback(hObject, eventdata, handles)
if handles.holdmods.Value
    X = evalin('base','Xnew');
else
    X = evalin('base','X');
end

if handles.shiftall.Value
    xs = str2double(handles.xshift.String);
    ys = str2double(handles.yshift.String);
    zs = str2double(handles.zshift.String);
    ts = str2double(handles.tshift.String);
    X = circshift(X,xs,1);
    X = circshift(X,ys,2);
    X = circshift(X,zs,3);
    X = circshift(X,ts,4);
else
    xs = str2double(handles.xshift.String);
    ys = str2double(handles.yshift.String);
    zs = str2double(handles.zshift.String);
    ts = str2double(handles.tshift.String);
    xp = str2double(handles.xpt.String);
    X(xp,:,:,:) = circshift(X(xp,:,:,:),zs,3);
end
assignin('base','Xnew',X);

% hObject    handle to shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function tshift_Callback(hObject, eventdata, handles)
% hObject    handle to tshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tshift as text
%        str2double(get(hObject,'String')) returns contents of tshift as a double


% --- Executes during object creation, after setting all properties.
function tshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zshift_Callback(hObject, eventdata, handles)
% hObject    handle to zshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zshift as text
%        str2double(get(hObject,'String')) returns contents of zshift as a double


% --- Executes during object creation, after setting all properties.
function zshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yshift_Callback(hObject, eventdata, handles)
% hObject    handle to yshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yshift as text
%        str2double(get(hObject,'String')) returns contents of yshift as a double


% --- Executes during object creation, after setting all properties.
function yshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xshift_Callback(hObject, eventdata, handles)
% hObject    handle to xshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xshift as text
%        str2double(get(hObject,'String')) returns contents of xshift as a double


% --- Executes during object creation, after setting all properties.
function xshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in shiftall.
function shiftall_Callback(hObject, eventdata, handles)
% hObject    handle to shiftall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of shiftall


% --- Executes on button press in align.
function align_Callback(hObject, eventdata, handles)
if handles.holdmods.Value
    X = evalin('base','Xnew');
else
X = evalin('base','X');
end
xs = str2double(handles.xshift.String);
if xs ~= 0
    xr = -xs:xs;
else
    xr = 0;
end
ys = str2double(handles.yshift.String);
if ys ~= 0
    yr = -ys:ys;
else
    yr = 0;
end
zs = str2double(handles.zshift.String);
ts = str2double(handles.tshift.String);
xp = str2double(handles.xpt.String);
yp = str2double(handles.ypt.String);
tp = str2double(handles.tpt.String);
Y = X;
for i = str2double(handles.tmin.String):str2double(handles.tmax.String)
%     for i = tp
%         figure;
        hold all;
    orig = squeeze(X(xp,yp,:,i));
    for j = xr
        for k = yr
            comp = squeeze(X(xp+j,yp+k,:,i));
%             plot(comp);
            test = corr(orig,comp);
           % disp(test);
            if abs(test) > 0.1  
                cor = xcorr(orig,comp);
%                 plot(cor);
                s = round(length(cor)/2)+find(cor == max(cor));
                if s~=0
%                     disp(s)
                    Y(xp+j,yp+k,:,i) = circshift(comp,s);
                end
            end
        end
    end
    multiWaitbar('Shifting Depths',i/str2double(handles.tmax.String));
end
assignin('base','Xnew',Y);
multiWaitbar('CLOSEALL');


% hObject    handle to align (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in plot5d.
function plot5d_Callback(hObject, eventdata, handles)

% H = handles.varselect.String{handles.varselect.Value};
M = handles.timeselect.String{handles.timeselect.Value};
x1 = str2double(handles.xmin.String);
x2 = str2double(handles.xmax.String);
xpt = str2double(handles.xpt.String);
y1 = str2double(handles.ymin.String);
y2 = str2double(handles.ymax.String);
ypt = str2double(handles.ypt.String);
z1 = str2double(handles.zmin.String);
z2 = str2double(handles.zmax.String);
zpt = str2double(handles.zpt.String);
tpt = str2double(handles.tpt.String);
tfpt = str2double(handles.tfpt.String);
tf1 = str2double(handles.tfmin.String);
tf2 = str2double(handles.tfmax.String);
if handles.plotmod.Value
    P = evalin('base','Xnew');
else
%     switch H
%         case 'Pressure'
%             P = evalin('base','pressure');
%         case 'Current'
%             P = evalin('base','I');
%         case 'AE'
%             P = evalin('base','AE');
%         case 'Lead'
%             P = evalin('base','L');
%         case 'Medium'
%             P = evalin('base','medium');
%         case 'Intensity'
%             P = evalin('base','sensor_data.p2');
%     end
 vname = handles.var_name.String;
    P = evalin('base',vname);
end
% if ndims(P) == 4
%     P = permute(P,[2 3 1 4 5]);
% end


% switch M
%     case 'Fast'
%         %                 P = P(x1:x2,y1:y2,z1:z2,tf1:tf2,tpt);
%         t = tpt;
%     case 'Slow'
%         P = permute(P,[1 2 3 5 4]);
%         %                 P = P(x1:x2,y1:y2,z1:z2,t1:t2,tfpt);
%         t = tpt;
%     case 'Depth'
%         P = permute(P,[1 2 4 3 5]);
%         %                 P = P(x1:x2,y1:y2,tf1:tf2,z1:z2,tpt);
%         t = zpt;
% end
c = [str2double(get(handles.cmin,'String')) str2double(get(handles.cmax,'String'))];
cmax = max(P(:))*c(2);
cmin = min(P(:))*c(1)*-1;
c = [cmin cmax];
hmap = get(handles.cmapmenu,'String');
h = hmap{get(handles.cmapmenu,'Value')};
if get(handles.cmapmenu,'Value') == 4
    clear h
    h = hotcoldDB;
elseif get(handles.cmapmenu,'Value') == 3
    clear h
    h = hotcold;
end
colormap(h);
% switch M
%     case {'Fast' 'Slow'}
axes(handles.axes1)
P2 = squeeze(P(x1:x2,y1:y2,zpt,tfpt)); % XY
if ~isnan(c)
    imagesc(P2,c)
else
    imagesc(P2)
end
xlabel('Y')
ylabel('X')
axes(handles.axes3)
P2 = squeeze(P(x1:x2,ypt,z1:z2,tfpt)); %XZ
if ~isnan(c)
    imagesc(P2,c)
else
    imagesc(P2)
end
xlabel('Z')
ylabel('X')
axes(handles.axes4)
P2 = squeeze(P(xpt,y1:y2,z1:z2,tfpt)); %YZ
if ~isnan(c)
    imagesc(P2,c)
else
    imagesc(P2)
end
xlabel('Z')
ylabel('Y')
axes(handles.axes7)
P2 = squeeze(P(x1:x2,ypt,zpt,tf1:tf2)); %XT
if ~isnan(c)
    imagesc(P2,c)
else
    imagesc(P2)
end
xlabel('T')
ylabel('X')
%     case 'Depth'
%         axes(handles.axes1)
%         P2 = squeeze(P(:,ypt,:,t)); % XT
%         axes(handles.axes3)
%         P2 = squeeze(P(xpt,:,:,t)); %YT
%         axes(handles.axes4)
%         P2 = squeeze(P(:,:,tfpt,t)); %XY
%         axes(handles.axes7)
%         P2 = squeeze(P(xpt,ypt,:,:))'; %ZT
% end

        
% hObject    handle to plot5d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in movie5d.
function movie5d_Callback(hObject, eventdata, handles)
if handles.plotmod.Value
    X = evalin('base','Xnew');
else
%  H = handles.varselect.String{handles.varselect.Value};
%     switch H
%         case 'Pressure'
%             X = evalin('base','pressure');
%         case 'Current'
%             X = evalin('base','I');
%         case 'AE'
%             X = evalin('base','AE');
%         case 'Lead'
%             X = evalin('base','L');
%         case 'Medium'
%             X = evalin('base','medium');
%         case 'Intensity'
%             X = evalin('base','sensor_data.p2');
%     end
 vname = handles.var_name.String;
    X = evalin('base',vname);
end

c = [str2double(get(handles.cmin,'String')) str2double(get(handles.cmax,'String'))];
cmax = max(X(:))*c(2);
cmin = min(X(:))*c(1)*-1;
c = [cmin cmax];
hmap = get(handles.cmapmenu,'String');
h = hmap{get(handles.cmapmenu,'Value')};
if get(handles.cmapmenu,'Value') == 4
    clear h
    h = hotcoldDB;
elseif get(handles.cmapmenu,'Value') == 3
    clear h
    h = hotcold;
end
colormap(h);

xpt = str2double(handles.xpt.String);
ypt = str2double(handles.ypt.String);
zpt = str2double(handles.zpt.String);
xr = str2double(handles.xmin.String):str2double(handles.xmax.String);
yr = str2double(handles.ymin.String):str2double(handles.ymax.String);
zr = str2double(handles.zmin.String):str2double(handles.zmax.String);
tr = str2double(handles.tfmin.String):str2double(handles.tfmax.String);
P1 = squeeze(X(xr,yr,zpt,tr));
P2 = squeeze(X(xr,ypt,zr,tr));
P3 = squeeze(X(xpt,yr,zr,tr));
for i = str2double(handles.tfmin.String):str2double(handles.tfmax.String)
    axes(handles.axes1)
    if ~isnan(c)
    imagesc(squeeze(P1(:,:,i)),c);
    else
         imagesc(squeeze(P1(:,:,i)));
    end
    xlabel('Depth');
    ylabel('Lateral');
    axes(handles.axes3)
    if ~isnan(c)
    imagesc(squeeze(P2(:,:,i))',c);
    else
          imagesc(squeeze(P2(:,:,i))');
    end
    xlabel('Depth');
    ylabel('Elevational');
    axes(handles.axes4)
    if ~isnan(c)
    imagesc(squeeze(P3(:,:,i))',c);
    else 
         imagesc(squeeze(P3(:,:,i))');
    end
    title(i);
    xlabel('Elevational');
    ylabel('Lateral');
end

% hObject    handle to movie5d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in varselect.
function varselect_Callback(hObject, eventdata, handles)
% hObject    handle to varselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns varselect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from varselect


% --- Executes during object creation, after setting all properties.
function varselect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to varselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in timeselect.
function timeselect_Callback(hObject, eventdata, handles)
% hObject    handle to timeselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns timeselect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from timeselect


% --- Executes during object creation, after setting all properties.
function timeselect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tfshift_Callback(hObject, eventdata, handles)
% hObject    handle to tfshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tfshift as text
%        str2double(get(hObject,'String')) returns contents of tfshift as a double


% --- Executes during object creation, after setting all properties.
function tfshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tfshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tfmin_Callback(hObject, eventdata, handles)
% hObject    handle to tfmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tfmin as text
%        str2double(get(hObject,'String')) returns contents of tfmin as a double


% --- Executes during object creation, after setting all properties.
function tfmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tfmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tfmax_Callback(hObject, eventdata, handles)
% hObject    handle to tfmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tfmax as text
%        str2double(get(hObject,'String')) returns contents of tfmax as a double


% --- Executes during object creation, after setting all properties.
function tfmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tfmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tfpt_Callback(hObject, eventdata, handles)
% hObject    handle to tfpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tfpt as text
%        str2double(get(hObject,'String')) returns contents of tfpt as a double


% --- Executes during object creation, after setting all properties.
function tfpt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tfpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function tfslide_Callback(hObject, eventdata, handles)
val = round(get(hObject,'Value'));
t = [str2double(handles.tfmin.String) str2double(handles.tfmax.String)];
if val > t(2)
    val = t(2)-1;
end
if val < t(1)
    val = t(1)+1;
end
handles.tfslide.Max = t(2);
handles.tfslide.Min = t(1);
set(handles.tfslide,'SliderStep',[1/t(2) 5/t(2)]);
set(handles.tfpt,'String',val);
set(handles.tfslide,'Value',val);
if handles.silder5d.Value
    plot5d_Callback(hObject, eventdata,handles);
else
    plotorig_Callback(hObject, eventdata, handles)
end
% hObject    handle to tfslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function tfslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tfslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in max5d.
function max5d_Callback(hObject, eventdata, handles)
% hObject    handle to max5d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 vname = handles.var_name.String;
    P = evalin('base',vname);
S = size(P);
%Sets min ranges to 1
set(handles.xmin,'String',1);
set(handles.ymin,'String',1);
set(handles.zmin,'String',1);
set(handles.tfmin,'String',1);
set(handles.tmin,'String',1);
set(handles.cmin,'String',num2str(min(P(:))));
set(handles.cmin,'String',-1);
set(handles.cmax,'String',1);
%Sets max ranges
set(handles.xmax,'String',S(1));
set(handles.ymax,'String',S(2));
switch length(S)
    case 2
        set(handles.zmax,'String',1);
        set(handles.tfmax,'String',1);
        set(handles.tmax,'String',1);
    case 3
        set(handles.zmax,'String',S(3));
        set(handles.tfmax,'String',1);
        set(handles.tmax,'String',1);
    case 4
        set(handles.zmax,'String',S(3));
        set(handles.tfmax,'String',S(4));
        set(handles.tmax,'String',1);
    case 5
        set(handles.zmax,'String',S(3));
        set(handles.tfmax,'String',S(4));
        set(handles.tmax,'String',S(5));
end


% --- Executes on button press in max_reg.
function max_reg_Callback(hObject, eventdata, handles)
if handles.plotmod.Value == 1
    X = evalin('base','Xnew');
else
    X = evalin('base','X');
end
S = size(X);
%Sets min ranges to 1
set(handles.xmin,'String',1);
set(handles.ymin,'String',1);
set(handles.zmin,'String',1);
set(handles.tmin,'String',1);
set(handles.cmin,'String',min(X(:)));
%Sets max ranges
set(handles.xmax,'String',S(1));
set(handles.ymax,'String',S(2));
set(handles.zmax,'String',S(3));
if length(S) >= 4
    set(handles.tmax,'String',S(4));
else
    set(handles.tmax,'String',1);
end
set(handles.cmax,'String',max(X(:)));
% hObject    handle to max_reg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in silder5d.
function silder5d_Callback(hObject, eventdata, handles)
% hObject    handle to silder5d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of silder5d



function plotlinex_Callback(hObject, eventdata, handles)
% hObject    handle to plotlinex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of plotlinex as text
%        str2double(get(hObject,'String')) returns contents of plotlinex as a double


% --- Executes during object creation, after setting all properties.
function plotlinex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotlinex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function plotliney_Callback(hObject, eventdata, handles)
% hObject    handle to plotliney (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of plotliney as text
%        str2double(get(hObject,'String')) returns contents of plotliney as a double


% --- Executes during object creation, after setting all properties.
function plotliney_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotliney (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plotline.
function plotline_Callback(hObject, eventdata, handles)
% hObject    handle to plotline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.plotmod.Value
    X = evalin('base','Xnew');
else
     vname = handles.var_name.String;
    X = evalin('base',vname);
end
%H = handles.varselect.String{handles.varselect.Value};
M = handles.timeselect.String{handles.timeselect.Value};
xpt = str2double(handles.xpt.String);
ypt = str2double(handles.ypt.String);
tfpt = str2double(handles.tfpt.String);
tpt = str2double(handles.tpt.String);
zpt = str2double(handles.zpt.String);

N = get(handles.plot2d,'Value');
if N == 1 %XZ
    p1 = str2double(get(handles.ypt,'String'));
    p2 = str2double(get(handles.tpt,'String'));
    x = str2double(get(handles.xmin,'String')):str2double(get(handles.xmax,'String'));
    y = str2double(get(handles.zmin,'String')):str2double(get(handles.zmax,'String'));
    imag = squeeze(X(x,p1,y,p2))';
elseif N == 2 %YZ
    p1 = str2double(get(handles.xpt,'String'));
    p2 = str2double(get(handles.tpt,'String'));
    x = str2double(get(handles.ymin,'String')):str2double(get(handles.ymax,'String'));
    y = str2double(get(handles.zmin,'String')):str2double(get(handles.zmax,'String'));
    imag = squeeze(X(p1,x,y,p2))';
elseif N == 3 %XY
    p1 = zpt;
    p2 = tfpt;
    p3 = tpt;
    x = str2double(get(handles.xmin,'String')):str2double(get(handles.xmax,'String'));
    y = str2double(get(handles.ymin,'String')):str2double(get(handles.ymax,'String'));
    if ndims(X) == 3
        imag = squeeze(X(x,y,p2)');
    elseif ndims(X) == 4
        imag = squeeze(X(x,y,p1,p2)');
    else
        imag = squeeze(X(x,y,p1,p2,p3));
    end

elseif N == 4 %TX
    p1 = str2double(get(handles.ypt,'String'));
    p2 = str2double(get(handles.zpt,'String'));
    x = str2double(get(handles.xmin,'String')):str2double(get(handles.xmax,'String'));
    y = str2double(get(handles.tmin,'String')):str2double(get(handles.tmax,'String'));
    imag = squeeze(X(x,p1,p2,y));
elseif N == 5 %TY
    p1 = str2double(get(handles.xpt,'String'));
    p2 = str2double(get(handles.zpt,'String'));
    x = str2double(get(handles.ymin,'String')):str2double(get(handles.ymax,'String'));
    y = str2double(get(handles.tmin,'String')):str2double(get(handles.tmax,'String'));
    imag = squeeze(X(p1,x,p2,y));
elseif N == 6 %TZ
    p1 = str2double(get(handles.xpt,'String'));
    p2 = str2double(get(handles.ypt,'String'));
    x = str2double(get(handles.zmin,'String')):str2double(get(handles.zmax,'String'));
    y = str2double(get(handles.tmin,'String')):str2double(get(handles.tmax,'String'));
    imag = squeeze(X(p1,p2,x,y));
end

if handles.op_env.Value
    imag = envelope(imag')';
end
if handles. op_invert.Value
    imag = imag*-1;
end
if handles.op_transpose.Value
    imag = imag';
end
axes(handles.axes1);
imagesc(imag);
hold on;
xcor = str2double(handles.plotlinex.String);
ycor = str2double(handles.plotliney.String);
plot(ones(size(imag,1))*xcor,'b--');
plot(ones(size(imag,1))*ycor','g--');




p = imag;
ymax = max(abs(p(:)));
yhalf = ymax/2;
if handles.auto_line.Value
    F = find(abs(p) == ymax,1);
    col = ceil(F/size(p,1));
    row = mod(F,size(p,1));
else
    col = xcor;
    row = ycor;
end
% axes(handles.axes1)
% imagesc(p);
pvert = p(:,col);
phorz = p(row,:);
pvert = resample(pvert,10,1);
axes(handles.axes3)
plot(pvert,'b','LineWidth',2);
if handles.op_fwhm_line.Value
    hold on;
    plot(ones(1,length(pvert))*yhalf,'r--')
    hold off
end
phorz = resample(phorz,10,1);
axes(handles.axes4)
hold off
cla reset
plot(phorz,'g','LineWidth',2);
if handles.op_fwhm_line.Value
    hold on
    plot(ones(1,length(phorz))*yhalf,'r--')
    hold off
end
%s2 = snum*10;
v_first = pvert(1:row*10);
v_second = pvert(row*10+1:end);
h_first = phorz(1:col*10);
h_second = phorz(col*10+1:end);
x1 = find(v_first <= yhalf,1,'last');
x2 = find(v_second <= yhalf,1)+row*10;
% kgrid = evalin('base','kgrid');
dy = str2double(handles.dy.String)/10;
dx = str2double(handles.dx.String)/10;
dz = str2double(handles.dz.String)/10;
dt = str2double(handles.dt.String)/10;
y1 = find(h_first <= yhalf,1,'last');
y2 = find(h_second <= yhalf,1)+col*10;
fwhm = (x2-x1)*dy;
fwhm2 = (y2-y1)*dx;
set(handles.peakx,'String',round(max(phorz),2));
set(handles.peaky,'String',max(pvert));
set(handles.fwhmy,'String',fwhm2);
set(handles.fwhmx,'String',fwhm);
set(handles.xloc,'String',col);
set(handles.yloc,'String',row);



% --- Executes on button press in loaddx.
function loaddx_Callback(hObject, eventdata, handles)
% hObject    handle to loaddx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
kgrid = evalin('base','kgrid');
set(handles.dx,'String',kgrid.dx*1e3)
set(handles.dy,'String',kgrid.dy*1e3)
set(handles.dz,'String',kgrid.dz*1e3)
set(handles.dt,'String',kgrid.dt*1e9);



function dx_Callback(hObject, eventdata, handles)
% hObject    handle to dx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dx as text
%        str2double(get(hObject,'String')) returns contents of dx as a double


% --- Executes during object creation, after setting all properties.
function dx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dy_Callback(hObject, eventdata, handles)
% hObject    handle to dy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dy as text
%        str2double(get(hObject,'String')) returns contents of dy as a double


% --- Executes during object creation, after setting all properties.
function dy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dz_Callback(hObject, eventdata, handles)
% hObject    handle to dz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dz as text
%        str2double(get(hObject,'String')) returns contents of dz as a double


% --- Executes during object creation, after setting all properties.
function dz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dt_Callback(hObject, eventdata, handles)
% hObject    handle to dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dt as text
%        str2double(get(hObject,'String')) returns contents of dt as a double


% --- Executes during object creation, after setting all properties.
function dt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in auto_line.
function auto_line_Callback(hObject, eventdata, handles)
% hObject    handle to auto_line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of auto_line


% --- Executes on button press in op_env.
function op_env_Callback(hObject, eventdata, handles)
set(handles.bband,'Value',0);
% hObject    handle to op_env (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of op_env


% --- Executes on button press in op_fwhm_line.
function op_fwhm_line_Callback(hObject, eventdata, handles)
% hObject    handle to op_fwhm_line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of op_fwhm_line


% --- Executes on button press in op_invert.
function op_invert_Callback(hObject, eventdata, handles)
% hObject    handle to op_invert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of op_invert


% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
axes(handles.axes1)
cla reset;
axes(handles.axes3)
cla reset;
axes(handles.axes4)
cla reset
axes(handles.axes7)
cla reset;
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function Plot_FFT(hObject,eventdata,handles)
a =1;
p = hObject.YData;
%p = p/max(abs(p));
axes(handles.axes7);
dx = str2double(handles.dx.String)/10;
Fs = 1/dx;
L = length(p);
x = (0:L-1)*dx;
if ismember('medium',evalin('base','who'))
    medium = evalin('base','medium');
    c = medium.sound_speed;
else
    c = 1500;
end
dt = str2double(handles.dt.String);
f = Fs*(0:(L/2))/L;
faxis = linspace(0,Fs,length(p));
P = fft(p);
P2 = abs(P/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
lim = round(length(p)/2);
%plot(f,P1);
fus = f*c/1000;
plot(fus,P1);
%plot(faxis(1:lim),abs(P(1:lim)));
xlim([0 Fs/10])

function Plot_FWHM(hObject,eventdata,handles)
pt = get(gca,'currentpoint');
x1 = round(pt(1,1));
y1 = round(pt(1,2));

p = hObject.CData;
pmean = mean(p(:));
p = p-pmean;
% if handles.op_env.Value
%     p = envelope(p')';
% end
% if handles.op_invert.Value
%     p = p*-1;
% end
if handles.dbc.Value
    p = real(20*log10(p./max(p(:))));
end
ymax1 = max(abs(p(:,x1)));
ymax2 = max(abs(p(y1,:)));
ymin2 = min(abs(p(y1,:)));
ymin1 = min(abs(p(:,x1)));
yhalf1 = (ymax1-ymin1)/2+ymin1;   
yhalf2 = (ymax2-ymin2)/2+ymin2;
col = x1;
row = y1;
% axes(handles.axes1)
% imagesc(p);
pvert = p(:,col);
phorz = p(row,:);
pvert = resample(pvert,10,1);
if handles.ext_fig.Value
    figure(26);
else
axes(handles.axes3)
end
plot(pvert,'b','LineWidth',2,'ButtonDownFcn',{@Plot_FFT,handles});
if handles.op_fwhm_line.Value
    hold on;
    plot(ones(1,length(pvert))*yhalf1,'r--')
    hold off
end
phorz = resample(phorz,10,1);
if handles.ext_fig.Value
    figure(27)
else
axes(handles.axes4)
end
hold off
cla reset
plot(phorz,'g','LineWidth',2,'ButtonDownFcn',{@Plot_FFT,handles});
if handles.op_fwhm_line.Value
    hold on
    plot(ones(1,length(phorz))*yhalf2,'r--')
    hold off
end
%s2 = snum*10;
v_first = pvert(1:row*10);
v_second = pvert(row*10+1:end);
h_first = phorz(1:col*10);
h_second = phorz(col*10+1:end);
x1 = find(v_first <= yhalf1,1,'last');
x2 = find(v_second <= yhalf1,1)+row*10;
% kgrid = evalin('base','kgrid');
dy = str2double(handles.dy.String)/10;
dx = str2double(handles.dx.String)/10;
dz = str2double(handles.dz.String)/10;
dt = str2double(handles.dt.String)/10;
y1 = find(h_first <= yhalf2,1,'last');
y2 = find(h_second <= yhalf2,1)+col*10;
fwhm = (x2-x1)*dy;
fwhm2 = (y2-y1)*dx;
set(handles.peakx,'String',round(max(phorz),2));
set(handles.peaky,'String',max(pvert));
set(handles.fwhmy,'String',fwhm2);
set(handles.fwhmx,'String',fwhm);
set(handles.xloc,'String',col);
set(handles.yloc,'String',row);

a = 2;


% --- Executes on button press in op_transpose.
function op_transpose_Callback(hObject, eventdata, handles)
% hObject    handle to op_transpose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of op_transpose


% --- Executes on button press in p_max.
function p_max_Callback(hObject, eventdata, handles)
% hObject    handle to p_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of p_max



function var_name_Callback(hObject, eventdata, handles)
% hObject    handle to var_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of var_name as text
%        str2double(get(hObject,'String')) returns contents of var_name as a double


% --- Executes during object creation, after setting all properties.
function var_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to var_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
vname = handles.var_name.String;
P = evalin('base',vname);
P2 = permute(P,[str2num(handles.p_order.String)]);
assignin('base',vname,P2);
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function p_order_Callback(hObject, eventdata, handles)
% hObject    handle to p_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p_order as text
%        str2double(get(hObject,'String')) returns contents of p_order as a double


% --- Executes during object creation, after setting all properties.
function p_order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in noiseselect.
function noiseselect_Callback(hObject, eventdata, handles)
% hObject    handle to noiseselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('Choose Signal');
[x2,y2] = ginput(1);
disp('Choose Noise');
[x,y] = ginput(1);
C = handles.axes1.Children.CData;

%Get Noise Parameters
xlims = floor(handles.axes1.XLim);
xmid = round(mean(xlims));
xrange = max(xlims)-min(xlims);
x4 = round(xrange/8);
xnoise = round(x)-x4:round(x)+x4;
ylims = floor(handles.axes1.YLim);
ymid = round(mean(ylims));
yrange = max(ylims)-min(ylims);
y4 = round(yrange/8);
ynoise = round(y)-y4:round(y)+y4;   
noise = mean(mean(abs(C(ynoise,xnoise))));

%Get Signal Parameters
x6 = round(xrange/12);
y6 = round(yrange/12);
xsig = round(x2)-x6:round(x2)+x6;   
ysig = round(y2)-y6:round(y2)+y6;   
sig = max(max(abs(C(ysig,xsig))));

%Calc and Disp SNR
snr = 20*log10(sig/noise);
set(handles.noise_x,'String',round(x))
set(handles.noise_y,'String',round(y))
set(handles.noise_snr,'String',round(snr,2))



function recon_P_Callback(hObject, eventdata, handles)
% hObject    handle to recon_P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of recon_P as text
%        str2double(get(hObject,'String')) returns contents of recon_P as a double


% --- Executes during object creation, after setting all properties.
function recon_P_CreateFcn(hObject, eventdata, handles)
% hObject    handle to recon_P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function recon_AE_Callback(hObject, eventdata, handles)
% hObject    handle to recon_AE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of recon_AE as text
%        str2double(get(hObject,'String')) returns contents of recon_AE as a double


% --- Executes during object creation, after setting all properties.
function recon_AE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to recon_AE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
Mn2 = handles.recon_AE.String;
P2 = handles.recon_P.String;
Mn = evalin('base',Mn2);

% H = hadamard(size(Mn,2));
%         V = squeeze(sum(V2,3));
 P = evalin('base',P2);
%V2 = zeros([size(P,1),size(P,2),size(P,3),min([size(P,4),size(Mn,1)]),size(Mn,4)]);
V2 = zeros([size(P,1),size(P,2),min([size(P,4),size(Mn,1)])]);
%Mn = permute(Mn,[3 1 2 4]);
V = zeros([size(P,1),size(P,2)]);
V3 = zeros([size(P,1),size(P,2),size(P,3)]);
% for k = 1:size(P,3)
%     P = evalin('base',P2);
%     Mn = evalin('base',Mn2);
%     for q = 1:size(P,3)
%         P(:,:,q,:) = P(:,:,q,:).*H(q,k);
% %         Mn = Mn.*H(k,:); 
%     end 
%     Mn = Mn.*H(k,:); %FULL HADAMARD
    for j = 1:size(P,3)
       for i  = 1:min([size(P,4),size(Mn,1)])
            % for k = 1:size(Mn,4)
            Temp(:,:,i) = P(:,:,j,i) .* Mn(i,j);
           % V2 = sum(V2,Temp);
      %      V2(:,:,j,i) = P(:,:,j,i) .* Mn(i,j); %Mn(i,j,1,k);
            % end

        end
        V2 = V2 + Temp;
        multiWaitbar('Reconstructing Image',j/size(P,3));
    end
    V3= V2; %SYNTHETIC
%     V3(:,:,k) = sum(V2,3); %FULL HADAMARD
%   %  V = V + V2;
%     multiWaitbar('Hadamard seqeuence',k/size(P,3));
%     multiWaitbar('Reconstructing Image',0);
% end
multiWaitbar('CLOSEALL');
V = sum(V3,3);
% assignin('base','V2',V2);
% V = squeeze(mean(V2,3));
% V = squeeze(mean(V,3))';
assignin('base','V2',V3); %SYNTHETIC = AxBxT;;; HADAMARD  = AxBxH
assignin('base','V',V);
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function reshape_box_Callback(hObject, eventdata, handles)
% hObject    handle to reshape_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of reshape_box as text
%        str2double(get(hObject,'String')) returns contents of reshape_box as a double


% --- Executes during object creation, after setting all properties.
function reshape_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reshape_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
P = handles.var_name.String;
P2 = evalin('base',P);
P3 = reshape(P2,[str2num(handles.reshape_box.String)]);
assignin('base',P,P3);
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton30.
function pushbutton30_Callback(hObject, eventdata, handles)
v = handles.con_menu.String{handles.con_menu.Value};
switch v
    case 'Hadamard'
        source = evalin('base','source');
        kgrid = evalin('base','kgrid');
        %  P = evalin('base',handles.recon_P.String);
        V = evalin('base',handles.had_output.String); %z,x,y,t
        if handles.tr_delays.Value
            delays = round(evalin('base','delays'));
            for i = 1:size(V,2)
            V(:,i) = circshift(V(:,i),delays(i));
            end
        end
        mask = source.p_mask;
        width = str2double(handles.con_dims.String);
        if ndims(mask) == 2

            y0 = ceil(find(mask,1)/kgrid.Nx);
            xline = mask(:,y0);
        elseif ndims(mask) == 3
       %     y0 = ceil(find(mask,1)/kgrid.Ny/kgrid.Nz);
         %   y0 = 
            y0 = 10;
            xline = mask(y0,:,round(size(mask,3)/2));
        end
        el_locs = find(xline);
        num_el = length(find(xline))/width;
        for i = 1:num_el
            for j = 1:width
                xloc(i,j) = el_locs((i-1)*width+j);
            end
        end
        xloc2 = xloc-kgrid.Nx/2;
        xloc3 = round(mean(xloc,2));
        if ndims(mask) == 2
        %x = (-kgrid.Nx/2:kgrid.Nx/2).*kgrid.dx*1000;
        x = kgrid.x_vec*1000;
        %     x = kgrid.x_vec./kgrid.dx;
        y = (kgrid.y_vec-min(kgrid.y_vec))*1000;
        %     y = (kgrid.y_vec-min(kgrid.y_vec))./kgrid.dy;
        elseif ndims(mask) == 3
                   x = kgrid.y_vec*1000;
        %     x = kgrid.x_vec./kgrid.dx;
        y = (kgrid.x_vec-min(kgrid.x_vec))*1000;
        end
        for k = 1:num_el %transmit
            for i = 1:length(x) %x points   
                if ndims(mask) == 2
                for j = 1:kgrid.Ny-y0 % y points
                    h(i,j,k) = abs(sqrt((x(i)-x(xloc3(k)))^2+(y(j))^2));
                    % h4(i,k) = abs(sqrt((x(i)-x(xloc3(k)).^2)));
                    %   h4(i,k) = abs(sqrt((x(i)-x(xloc3(k)))^2));
                    co(i,j,k) = abs(y(j)/h(i,j,k));
                end
                elseif ndims(mask) == 3
                      for j = 1:kgrid.Nx-y0
                          h(i,j,k) = abs(sqrt((x(i)-x(xloc3(k)))^2+(y(j))^2));
                          co(i,j,k) = abs(y(j)/h(i,j,k));
                      end
                end
            end
        end
        co(isnan(co)) = 0;
        medium = evalin('base','medium');
        dep = kgrid.Ny*kgrid.dy*1000;
       % time = kgrid.dt*medium.sound_speed*kgrid.Nt*1000; %in mm
       
        time = kgrid.dt*medium.sound_speed(1)*size(V,1)*1000; %in mm
        T = linspace(0,time,size(V,1));
        ae_size = [size(h) size(V,3) size(V,4)];
        ae = zeros(ae_size);
        h2 = round(h./(mean([kgrid.dy kgrid.dx])*1000));
        t2 = round(T./(mean([kgrid.dy kgrid.dx])*1000));
        tic
        for p = size(V,4)
            for n = 1:size(V,3)
                for k = 1:num_el
                    for t = 1:size(V,1)
                        %disp(t);
                        temp = find(h2(:,:,k) == t2(t));
                        if length(temp) > 0
                            for i = 1:length(temp)
                                locy(i) = ceil(temp(i)/ae_size(1));
                                locx(i) = mod(temp(i),ae_size(1));
                                if locx(i) == 0
                                    locx(i) = ae_size(1);
                                end
                             %   ae(locx(i),locy(i),k) = V(t,k);
                               ae(locx(i),locy(i),k) = V(t,k).*co(locx(i),locy(i),k);
                            end
                        end
                    end
                    multiWaitbar('Constructing Image',k/num_el);
                end
                multiWaitbar('Elevational',n/size(V,3));
            end
            multiWaitbar('Over time',p/size(V,4));
        end
        toc
            multiWaitbar('CLOSEALL');
        assignin('base','V2',ae)
        V = sum(ae,3);
        assignin('base','V',V);
    case 'Planes'
        %%%%%%%%%%%%%%%%%% OLD CODE START
   source = evalin('base','source');
        kgrid = evalin('base','kgrid');
        %  P = evalin('base',handles.recon_P.String);
        V = evalin('base',handles.had_input.String);
        mask = source.p_mask;
        t_all = str2num(handles.con_dims.String);
        theta = t_all(1);
        dtheta = t_all(2);
        y0 = ceil(find(mask,1)/kgrid.Nx);
        xline = mask(:,y0);
        el_locs = find(xline);
        num_el = length(find(xline));
        H = str2num(handles.had_output.String);
        els = H(1);
        width = H(2);
        if width > 1
            if mod(width,2) ~= 0
                width = width-1;
            end
            for i = 1:els
                x_loc(i) = el_locs(1+width/2+width*(i-1));
            end
        else
            x_loc = el_locs;
        end
        xloc2 = (x_loc-kgrid.Nx/2).*kgrid.dx.*1000;
        xloc3 = round(x_loc);
        %x = (-kgrid.Nx/2:kgrid.Nx/2).*kgrid.dx*1000;
        x = kgrid.x_vec*1000;
        %     x = kgrid.x_vec./kgrid.dx;
        y = (kgrid.y_vec-min(kgrid.y_vec))*1000;
        n_fire = theta/dtheta+1;
        angles = (0:dtheta:theta) - (theta/2);
        %%%%%%%%%%%%% Calc Delays
        N = theta/dtheta+1;
        t = (0:dtheta:theta) - (theta/2);
        for i = 1:N
            %         for j = 1:length(xloc2)
            %             %  if t(i) <= 0
            %             %                 delays(i,j) = tand(t(i)).*abs(xloc2(1)-xloc2(j));
            %             delays(i,j) = tand(t(i)).*xloc2(end - j+1);
            %             % else
            %             %    delays(i,j) = tand(t(i)).*abs(xloc3(end)-xloc3(j));
            %             % end
            %
            %         end
            delays(i,:) = tand(t(i)).*flipud(xloc2); %mm
        end
        delays = delays-min(delays(:));
        delays  = fliplr(delays);
%         for i = 1:N
%             delays(i,:) = delays(i,:)-min(delays(i,:));
%         end
%         if mod(num_el,2) == 0
%             a = num_el/2;
%             b = a+1;
%             c = delays(:,[a b]);
%             d = sum(c,2)/2;
%             e = max(d);
%             f = abs(d-e);
%            delays = delays+f;
%         elseif mod(num_el,2) ~= 0
%             a = ceil(num_el/2);
%             d = delays(:,a);
%             e = max(d);
%             f = abs(d-e);
%             delays = delays+f;
%         end

        %%%%%%%%%%%%%%%%%%%%% burst offset method %%%%%%%%%%%%%%%%%%%%
        %     t0 = 1;
        %
        %      burst_offset = t0 + element_space * sin(steering_angle*pi/180)/...
        %                                 (ss * kgrid.dt);
        %                             burst_offset = burst_offset - min(burst_offset);
        %%%%%%%%%%%%%%%%%%%%%%%%ADJUSTED FROM FOCUS%%%%%%%%%%%%%%%%%%%%
%         medium = evalin('base','medium');
%          xstep = kgrid.dx*1000;%param.velmex.XDist/(param.velmex.XNStep-1);
%             ystep = kgrid.dy*1000; %ax.depth(end)-ax.depth(end-1);
%             tstep = kgrid.dt*1000; %ms
%             zstep = medium.sound_speed(1)*tstep; %These are all in mm
% 
% 
%          %   z = f(1).*kgrid.dy*1000; %mm
%           %  x = kgrid.x(xloc).*1000;
%             Y = zeros(size(V));
%             Y2 = zeros(kgrid.Nx,kgrid.Ny-y0)';
%             Y2 = zeros(kgrid.Nx,kgrid.Ny);
%             ng = zeros(size(V));
%             xmid = median(x);
% 
% 
% 
%             depth = 0:zstep:zstep*size(V,1);
%             depth = depth(1:end-1); %mm
% 
%         for q = 1:length(angles)
% 
%             for i = 1:size(V,1)
%                 for j = 1:els
%                     newx(i,j) = xloc2(j)+(depth(i)-delays(q,j))*sind(angles(q)); %mm
%                     newz(i,j) = (depth(i)-delays(q,j))*cosd(angles(q));
% %                        newx(i,j) = xloc2(j)+(depth(i))*sind(angles(q)); %mm
% %                     newz(i,j) = (depth(i))*cosd(angles(q));
%                 end
%             end
%             newx2 = round(newx./xstep); %samp
%             newz2 = round(newz./zstep);
%             new_depth = round(depth./zstep);
%             % D = Y;
%             % for i = 1:length(dx)
%             %     for j = 1:length(depth)
%             %
%             %         xx = dx(i);
%             %         x2 = find(dx == xx);
%             %         xx2 = find(dx == newx2(j,i));
%             %         zz2 = find(new_depth == newz2(j,i));
%             %         if Y(zz2,xx2) > 0
%             %             Y(zz2,xx2) = (Y(zz2,xx2) + V(j,x2))/2;
%             %         else
%             %             Y(zz2,xx2) = V(j,x2);
%             %         end
%             %
%             %     end
%             %     multiWaitbar('Converting',i/length(dx));
%             % end
%             %Y2 = Y2';
%             for i = 1:els
%                 for j = 1:length(depth)
% 
%                     xx = xloc2(i);  %mm loc
%                     x2 = find(xloc2 == xx); %transmit loc
%                     xx2 = find(round(kgrid.x_vec./kgrid.dx) == newx2(j,i));
%                     zz2 = find(round((kgrid.y_vec-min(kgrid.y_vec))./(zstep/1000)) == newz2(j,i));
% 
%                     if Y2(xx2,zz2) > 0
%                         Y2(xx2,zz2) = (Y2(xx2,zz2) + V(j,q))/2;
%                     else
%                         Y2(xx2,zz2) = V(j,q);
%                     end
%                 end
%                 multiWaitbar('Converting',i/els);
%             end
%             multiWaitbar('Compiling the full transducer',q/length(angles));
%             Y3(:,:,q) = Y2;
%             Y2 = zeros(kgrid.Nx,kgrid.Ny);
%         end
%         multiWaitbar('CLOSEALL');
%      %   Y3 = sum(Y3,3);
% assignin('base','Z',Y3);




        %%%%%%%%%%%%%%%%%%%%%%% OLD CODE FROM HADAMARD HERE %%%%%%%
        for m = 1:n_fire
            %delays = get_plane_delays(xloc3,theta,dtheta);
            for k = 1:els
                for i = 1:length(x)
                    for j = 1:kgrid.Ny-y0
                        if abs(sqrt((x(i)-x(xloc3(k)))^2+(y(j))^2)) > delays(m,k)
                            h(i,j,k,m) = abs(sqrt((x(i)-x(xloc3(k)))^2+(y(j))^2) - delays(m,k));
%                             if m == 1
%                                 w(i,j,k) = y(j)/abs(sqrt((x(i)-x(xloc3(k)))^2+(y(j))^2));%weight
%                             end
                        else
                            h(i,j,k,m) = 0;
%                             if m == 1
%                                 w(i,j,k) = 0;
%                             end
                        end
                        if m == 1
                            w(i,j,k) = y(j)/abs(sqrt((x(i)-x(xloc3(k)))^2+(y(j))^2));%weight
                        end
                        %                         if angles(m) < 0
                        %                             if delays(m,k) >= abs(sqrt((x(i)-x(xloc3(k)))^2+(y(j))^2))
                        %                                 h(i,j,k,m) = 0;
                        %                             else
%                                 h(i,j,k,m) = abs(sqrt((x(i)-x(xloc3(k)))^2+(y(j))^2)) - delays(m,k);
%                             end   
%                         elseif angles(m) == 0
%                             h(i,j,k,m) = abs(sqrt((x(i)-x(xloc3(k)))^2+(y(j))^2));
%                         elseif angles(m)> 0
%                             if delays(m,k) >= abs(sqrt((x(i)-x(xloc3(k)))^2+(y(j))^2))
%                                 h(i,j,k,m) = 0;
%                             else
%                                 h(i,j,k,m) = abs(sqrt((x(i)-x(xloc3(k)))^2+(y(j))^2)) - delays(m,k);
%                             end
%                         end
                    end
                end

            end
            multiWaitbar('Calculating Plane Angles',m/n_fire);
        end
        w(isnan(w)) = 0;
        for i = 1:size(delays,1)
            if max(delays(i,:)) > 0
                h(:,:,:,i) = h(:,:,:,i) + max(delays(i,:));
            end
        end

        medium = evalin('base','medium');
        dep = kgrid.Ny*kgrid.dy*1000;
        time = kgrid.dt*medium.sound_speed*kgrid.Nt*1000; %in mm
        T = linspace(0,time,size(V,1));
        ae_size = size(h);
        ae = zeros(ae_size);
        h2 = round(h./(mean([kgrid.dy kgrid.dx])*1000));
        t2 = round(T./(mean([kgrid.dy kgrid.dx])*1000));
        %    h3 = squeeze(mean(h2,3));
        %    ae_size = size(h3);
        %    ae = zeros(ae_size);

     %   for m = 1:n_fire
         % for k = 1:num_el
%             for t = 1:size(V,1)
%                 %disp(t);
%                 temp = find(h2 == t2(t));
%                 % z = floor(temp./(ae_size(1)*ae_size(2)));
%                 if length(temp) > 0
%                     for i = 1:length(temp)
%                         %I can introduce locz as the element
%                         [locy(i), locx(i), locz(i), locm(i)] = ind2sub(size(ae),temp(i));
%                         %                         locy(i) = ind2sub() %ceil(temp(i)/ae_size(1));
%                         %                         locx(i) = mod(temp(i),ae_size(1));
%                         %                         locz(i) =
%                         %                         locm(i) =
%                         %                         if locx(i) == 0
%                         %                             locx(i) = ae_size(1);
%                         %                         end
%                         ae(locx(i),locy(i),locz(i),locm(i)) = V(t,locm(i));
%                     end
%                 end
%                 multiWaitbar('Constructing Image',t/size(V,1));
% 
%                 %             multiWaitbar('Constructing Image Transmit',m/n_fire);
%             end


for m = 1:n_fire
    for k = 1:els
        for t = 1:size(V,1)
            %disp(t);
            temp = find(h2(:,:,k,m) == t2(t));
%             [w1 w2] = ind2sub([length(x) kgrid.Ny-y0],temp);
          %  weight = w(w1,w2,k,m);

            % z = floor(temp./(ae_size(1)*ae_size(2)));
            if length(temp) > 0
                for i = 1:length(temp)
                    %I can introduce locz as the element

                    locy(i) = ceil(temp(i)/ae_size(1));
                    locx(i) = mod(temp(i),ae_size(1));
                    if locx(i) == 0
                        locx(i) = ae_size(1);
                    end
%                     ae(locx(i),locy(i),k,m) = V(t,m);
                   ae(locx(i),locy(i),k,m) = V(t,m).*w(locx(i),locy(i),k);
                end
            end
        end
        multiWaitbar('Constructing Image Element',k/els);
    end
    multiWaitbar('Constructing Image Transmit',m/n_fire);
end

        multiWaitbar('CLOSEALL');
        assignin('base','V2',ae)
        V = squeeze(sum(ae,3));
        V = mean(V,3);
        assignin('base','V',V);




        %%%%%%%%%%%%%%%%% OLD CODE END
%         source = evalin('base','source');
%         kgrid = evalin('base','kgrid');
%         %  P = evalin('base',handles.recon_P.String);
%         V = evalin('base',handles.had_output.String);
%         mask = source.p_mask;
%         t_all = str2num(handles.con_dims.String);
%         theta = t_all(1);
%         dtheta = t_all(2);
%         y0 = ceil(find(mask,1)/kgrid.Nx);
%         xline = mask(:,y0);
%         el_locs = find(xline);
%         num_el = length(find(xline));
%         xloc2 = (el_locs-kgrid.Nx/2).*kgrid.dx.*1000;
%         xloc3 = round(mean(el_locs,2));
%         %x = (-kgrid.Nx/2:kgrid.Nx/2).*kgrid.dx*1000;
%         x = kgrid.x_vec*1000;
%         %     x = kgrid.x_vec./kgrid.dx;
%         y = (kgrid.y_vec-min(kgrid.y_vec))*1000;
%         n_fire = theta/dtheta+1;
%         angles = (0:dtheta:theta) - (theta/2);
%         %%%%%%%%%%%%% Calc Delays
%         N = theta/dtheta+1;
%         t = (0:dtheta:theta) - (theta/2);
%         xone = find(xline,1);
%         xtwo = length(xline) - find(flip(xline),1)+1;
%         xline2 = xone:xtwo;
%         new_el = length(xline2);
%         for i = 1:N
%             %         for j = 1:length(xloc2)
%             %             %  if t(i) <= 0
%             %             %                 delays(i,j) = tand(t(i)).*abs(xloc2(1)-xloc2(j));
%             %             delays(i,j) = tand(t(i)).*xloc2(end - j+1);
%             %             % else
%             %             %    delays(i,j) = tand(t(i)).*abs(xloc3(end)-xloc3(j));
%             %             % end
%             %
%             %         end
%             delays(i,:) = tand(t(i)).*flipud(xloc2); % delays in mm
%         end
%         for i = 1:N
%             delays(i,:) = delays(i,:)-min(delays(i,:));
%         end
%      %   delays = delays';
%     for i = 1:N
%         delays2(i,:) = linspace(delays(i,1),delays(i,end),new_el);
%     end
% delays = delays2;
%         %%%%%%%%%%%%%%%%%%%%% burst offset method %%%%%%%%%%%%%%%%%%%%
%         %     t0 = 1;
%         %
%         %      burst_offset = t0 + element_space * sin(steering_angle*pi/180)/...
%         %                                 (ss * kgrid.dt);
%         %                             burst_offset = burst_offset - min(burst_offset);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %delays = fliplr(delays);
%         %%%%%%%%%%%%%
%         h = zeros(length(x),kgrid.Ny-y0,new_el,n_fire);
%         for m = 1:n_fire
%             %delays = get_plane_delays(xloc3,theta,dtheta);
%           %  for k = 1:new_el %changed from num_el
%                 %                 for i = 1:length(x)
%                 %                     for j = 1:kgrid.Ny-y0
%               %  if k == 1 || k == new_el
%                    % if angles(m) ~= 0
%                         for i = 1:length(x)
%                             for j = 1:kgrid.Ny-y0
%                                 %if delays(m,k) >= abs(sqrt((x(i)-x(xline2(k)))^2+(y(j))^2)) %changed all xloc3 to xline2
%                                %     h(i,j,k,m) = 0;
%                               %  else
%                                     h(i,j,:,m) = (abs(sqrt((x(i)-x(xline2(:))).^2+(y(j))^2)))' - delays(m,:);
%                                % end
%                             end
%                         end
%                         multiWaitbar('Calculating Plane Angles',m/n_fire);
%         end
% %                     elseif angles(m) == 0
% %                         for i = 1:length(x)
% %                             for j = 1:kgrid.Ny-y0
% %                                 h(i,j,k,m) = abs(sqrt((x(i)-x(xline2(k)))^2+(y(j))^2));
% %                             end
% %                         end
% %                         %                             elseif angles(m)> 0
% %                         %                                 if delays(m,k) >= abs(sqrt((x(i)-x(xline2(k)))^2+(y(j))^2))
% %                         %                                     h(i,j,k,m) = 0;
% %                         %                                 else
% %                         %                                     h(i,j,k,cm) = abs(sqrt((x(i)-x(xline2(k)))^2+(y(j))^2)) - delays(m,k);
% %                         %                                 end
% %                         %                             end
% %                     else
% %                         for i = 1:length(x)
% %                             if i == xline2(k)
% %                                 for j = 1:kgrid.Ny-y0
% %                                     h(i,j,k,m) = abs(sqrt((x(i)-x(xline2(k)))^2+(y(j))^2))-delay(m,k);
% %                                 end
% %                             end
% %                         end
% %                     end
% %                 end
% %                 multiWaitbar('Calculating Soemthing',k/new_el);
% %             end
% %             multiWaitbar('Calculating Plane Angles',m/n_fire);
% %         end
%         for i = 1:size(delays,1)
%             if max(delays(i,:)) > 0
%                 h(:,:,:,i) = h(:,:,:,i) + max(delays(i,:));
%             end
%         end
% 
%         medium = evalin('base','medium');
%         dep = kgrid.Ny*kgrid.dy*1000;
%         time = kgrid.dt*medium.sound_speed*size(V,1)*1000; %in mm
%         T = linspace(0,time,size(V,1));
% %         ae_size = size(h);
% %         ae = zeros(ae_size);
%         h2 = round(h./(mean([kgrid.dy kgrid.dx])*1000));
%         t2 = round(T./(mean([kgrid.dy kgrid.dx])*1000));
%         %    h3 = squeeze(mean(h2,3));
%         %    ae_size = size(h3);
%         %    ae = zeros(ae_size);
% h3 = mean(h2,4);
% h_loc = round(linspace(1,new_el,n_fire-1));
% h4 = h3(:,:,h_loc);
%      ae_size = size(h4);
%         ae = zeros(ae_size);
% 
% 
% for k = 1:num_el
%     for t = 1:size(V,1)
%         %disp(t);
%         temp = find(h4(:,:,k) == t2(t));
%         if length(temp) > 0
%             for i = 1:length(temp)
%                 locy(i) = ceil(temp(i)/ae_size(1));
%                 locx(i) = mod(temp(i),ae_size(1));
%                 if locx(i) == 0
%                     locx(i) = ae_size(1);
%                 end
%                 ae(locx(i),locy(i),k) = V(t,k);
%             end
%         end
%     end
%     multiWaitbar('Constructing Image',k/num_el);
% end
% multiWaitbar('CLOSEALL');
% assignin('base','V2',ae)
% V = sum(ae,3);
% assignin('base','V',V);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         for m = 1:n_fire
%             for k = 1:num_el
%                 for t = 1:size(V,1)
%                     %disp(t);
%                     temp = find(h2(:,:,k,m) == t2(t));
%                     % z = floor(temp./(ae_size(1)*ae_size(2)));
%                     if length(temp) > 0
%                         for i = 1:length(temp)
%                             %I can introduce locz as the element
% 
%                             locy(i) = ceil(temp(i)/ae_size(1));
%                             locx(i) = mod(temp(i),ae_size(1));
%                             if locx(i) == 0
%                                 locx(i) = ae_size(1);
%                             end
%                             ae(locx(i),locy(i),k,m) = V(t,m);
%                         end
%                     end
%                 end
%                 multiWaitbar('Constructing Image Element',k/num_el);
%             end
%             multiWaitbar('Constructing Image Transmit',m/n_fire);
%         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    for m = 1:n_fire
        %        for t = 1:size(V,1)
        %            %disp(t);
        %            temp = find(h2(:,:,:,m) == t2(t));
        % %            z = floor(temp./(ae_size(1)*ae_size(2)));
        % %            temp = temp-z*(ae_size(1)*ae_size(2));
        %            if length(temp) > 0
        %                for i = 1:length(temp)
        %                    %I can introduce locz as the element
        %                    locz(i) = ceil(temp(i)./(ae_size(1)*ae_size(2)));
        %                    locy(i) = ceil(temp(i)/ae_size(1)./locz(i));%-z(i);
        %                    locx(i) = mod(temp(i),ae_size(1));%-z(i);
        %                    if locx(i) == 0
        %                        locx(i) = ae_size(1);
        %                    end
        %                    ae(locx(i),locy(i),locz(i),m) = V(t,m);
        %
        %                end
        %            end
        %             multiWaitbar('Pooping',t/size(V,1));
        %        end
        %        multiWaitbar('Constructing Image Transmit',m/n_fire);
        %    end
% 
% 
%         multiWaitbar('CLOSEALL');
%         assignin('base','V2',ae)
%         V = squeeze(sum(ae,3));
%         assignin('base','V',V);

    case 'Focus'
         source = evalin('base','source');
        kgrid = evalin('base','kgrid');
        V = evalin('base',handles.had_input.String);
        ad = str2double(handles.had_output.String); %added delay
%         mask = source.p_mask;
%       %  width = str2double(handles.con_dims.String);
%         y0 = ceil(find(mask,1)/kgrid.Nx);
%         xline = mask(:,y0);
%         el_locs = find(xline);
%       %  num_el = length(find(xline))/width;
        f = str2num(handles.con_dims.String);
%         dx = f(2);
%         lat = f(1);
%         n = size(V,2);
%         foci = (-(dx/2*(n-1)):dx:(dx/2*(n-1)))+lat;
%         depth = f(3);
%         xloc2 = el_locs-kgrid.Nx/2;
%         xloc3 = round(mean(el_locs,2));
%         %x = (-kgrid.Nx/2:kgrid.Nx/2).*kgrid.dx*1000;
%         x = kgrid.x_vec*1000;
%         %     x = kgrid.x_vec./kgrid.dx;
%         y = (kgrid.y_vec-min(kgrid.y_vec))*1000;
%         %     y = (kgrid.y_vec-min(kgrid.y_vec))./kgrid.dy;
%         xmid = round(median(xloc3));


        medium = evalin('base','medium');
xstep = kgrid.dx*1000;%param.velmex.XDist/(param.velmex.XNStep-1);
ystep = kgrid.dy*1000; %ax.depth(end)-ax.depth(end-1);
tstep = kgrid.dt*1000; %ms
zstep = medium.sound_speed(1)*tstep; %These are all in mm

mask = source.p_mask;
y0 = ceil(find(mask,1)/kgrid.Nx);
xline = mask(:,y0);
el_locs = find(xline);


z = f(1).*kgrid.dy*1000; %mm
x = kgrid.x(el_locs).*1000;
Y = zeros(size(V));
Y2 = zeros(kgrid.Nx,kgrid.Ny-y0)';
Y2 = zeros(kgrid.Nx,kgrid.Ny);
ng = zeros(size(V));
xmid = median(x);
%fpoints = 
   d_focus = z - y0*kgrid.dy*1000; %focus for delay calc

for q = 1:length(el_locs)

    xmid(q) = x(q);

    if mod(size(V,2),2) ~= 0
        Nx = (size(V,2)-1)/2;
    else
        Nx = size(V,2)/2;
    end
    fpoints  = (-Nx:1:Nx).*kgrid.dx.*1000.*f(2);
    dx(:,q) = (-Nx:1:Nx).*kgrid.dx.*1000.*f(2)-xmid(q); %mm Gives how far in x the focal point is relative to the chosen element

end


    for i = 1:length(x) % Each element
        for j = 1:length(fpoints) % Each focal point
            r(j,i) = sqrt(dx(j,i)^2+d_focus^2); % Element i and focal points j
          %   r(i,j,q) = sqrt((dx-xmid+x(i)).^2+(d_focus).^2);
        end
    end

%     r2 = abs(max(r,[],1)-r);
     r2 = abs(max(r(:))-r); %Delay added by focal points
    theta = atand(dx./d_focus);
    depth = 0:zstep:zstep*size(V,1);
    depth = depth(1:end-1); %mm This is actually radius

    %     delays = sqrt()

    for q = 1:length(el_locs)
        for i = 1:size(V,1)
            for j = 1:length(fpoints)
                newx(i,j,q) = xmid(q)+(depth(i)-r2(j,q))*sind(theta(j,q)); %mm
                newz(i,j,q) = (depth(i)-r2(j,q)).*cosd(theta(j,q))-ad;
            end
        end
    end
    newx2 = round(newx./xstep); %samp
    newz2 = round(newz./zstep);
    new_depth = round(depth./zstep);
    % D = Y;
    % for i = 1:length(dx)
    %     for j = 1:length(depth)
    %
    %         xx = dx(i);
    %         x2 = find(dx == xx);
    %         xx2 = find(dx == newx2(j,i));
    %         zz2 = find(new_depth == newz2(j,i));
    %         if Y(zz2,xx2) > 0
    %             Y(zz2,xx2) = (Y(zz2,xx2) + V(j,x2))/2;
    %         else
    %             Y(zz2,xx2) = V(j,x2);
    %         end
    %
    %     end
    %     multiWaitbar('Converting',i/length(dx));
    % end
    %Y2 = Y2';
    for q = 1:length(el_locs)
        for i = 1:length(dx)
            for j = 1:length(depth)
                %
                %             xx = dx(i);  %mm loc
                %             x2 = find(dx == xx); %transmit loc
                                xx2 = find(round(kgrid.x_vec./kgrid.dx) == newx2(j,i,q));
                                zz2 = find(round((kgrid.y_vec-min(kgrid.y_vec))./(zstep/1000)) == newz2(j,i,q));
%                 xx2 = newx2(j,i,q);
%                 zz2 = newz2(j,i,q);
                if zz2 < 1
                    zz2 = 1;
                elseif zz2 > 300
                    zz2 = 300;
                end

                if Y2(xx2,zz2) > 0
                    Y2(xx2,zz2) = (Y2(xx2,zz2) + V(j,i))/2;
                else
                    Y2(xx2,zz2) = V(j,i);
                end
            end
            multiWaitbar('Converting',i/length(dx));
        end
        multiWaitbar('Compiling the full transducer',q/length(x));
        Y3(:,:,q) = Y2;
        Y2 = zeros(kgrid.Nx,kgrid.Ny);
    end

multiWaitbar('CLOSEALL');
%assignin('base','Y3',Y)
assignin('base','Y2',Y3);
Y = sum(Y3,3);
assignin('base','Y',Y);






end


%     %%%%%%%%%%%
%     new_dims = str2num(handles.con_dims.String);
%     ScanPt = 1;
%     [file path] = uigetfile('*.dat');
%     fname = fullfile(path,file);
% 
%     froot = fname(1:end-9);
%     fparam = [froot '_info.mat'];
%     fLF = [froot '_LF_Avg.dat'];
%     fHF = [froot '_HF_Avg.dat'];
%     m = ScanPt; %Do I need this?
% 
%     %  HF = read_hfdata(fHF,m,param,p);
%     hffile = fHF;
%     ScanPt= m;
% 
%     blk_idx = ScanPt;
% 
%     %     if exist(hffile,'file')
%     fid = fopen(fname,'rb');
%     blk_idx = 1;
%     %     else
%     %         hf_avg_file = regexprep(hffile,'_P[0-9]{4,4}_','_');
%     %         fid = fopen(hf_avg_file,'rb');
%     %     end
% 
%     if fid < 0
%         error('Cannot open HF data file');
%     end
% 
%     sztype = fgets(fid,fread(fid,1,'int32'));
%     if(strcmpi(sztype(1:5),'ucsdi') == 0)
%         fclose(fid);
%         return;
%     end
%     n = fread(fid,1,'int32');
%     % dsize = fread(fid,[1,n],'int32');
%     dsize = str2num(handles.con_dims.String);
%     blk_size = prod(dsize)*4;
%     fseek(fid,blk_size*(ScanPt-1),'cof');
%     data = fread(fid,[prod(dsize),1],'single');
%     % data = fread(fid,[numel(dsize),1],'single');
% 
%     assignin('base','data',data);
% 
% 
%     fclose(fid);



    
% function [data,parm] = read_hfdata_v161108(fid,ScanPt,param,p)
%
%     % nPos = ftell(fid);
%     n = fread(fid,1,'int32');
%     dsize = fread(fid,[1,n],'int32');
%     if numel(dsize)<4
%         if param.post.ind || param.sweep
%             %dsize(4) = param.Scan.Avg;
%             dsize(4) = 1;
%         end
%     end
%     % blk_size = (length(dsize) + 1)*4 + prod(dsize)*4;
%     blk_size = prod(dsize)*4;
%     fseek(fid,blk_size*(ScanPt-1),'cof');
%     data = fread(fid,[prod(dsize),1],'single');
%     if param.post.ind || param.sweep
%         tsize = dsize(1)*dsize(2)*dsize(3)*dsize(4); %Slow Time, Channels, Fast Time, Averages
%     else
%         tsize = dsize(1)*dsize(2)*dsize(3);
%     end
%     if tsize ~= length(data)
%         data = zeros(tsize,1);
%         disp(['Missing point ' num2str(ScanPt)]);
%     end
%     if param.post.ind
%         data = reshape(data,[dsize(3),dsize(2),dsize(1),dsize(4)]);
%         data = permute(data,[3,1,4,2]);
%         data = data(:,:,:,p);
%     elseif param.sweep
%         d2 = reshape(data,[round(param.Daq.HF.Samples),round(param.Scan.Sweep_Tpt*param.Scan.Duration_ms)*param.Scan.Xpt*param.Scan.Ypt,1]); %param.Scan.Avg]);
%         d3 = mean(d2,3);
%         %     for i = 1:param.Scan.Xpt
%         for j = 1:round(param.Scan.Duration_ms*param.Scan.Sweep_Tpt)
%             HF(:,1:param.Scan.Xpt,j) = d3(:,(param.Scan.Xpt*(j-1)+1:param.Scan.Xpt*j));
%         end
%         %         for i = 1:param.Scan.Xpt
%         %             HF(:,i,1:param.Scan.Sweep_Tpt) = d3(:,(param.Scan.Sweep_Tpt*(i-1)+1):param.Scan.Sweep_Tpt*i);
%         %         end
%         % for i = 1:param.Scan.Xpt
%         %         for j = 1:param.Scan.Sweep_Tpt
%         %             HF(i,:,j,:) = d2(:,i+(param.Scan.Xpt*(j-1)),:);
%         %         end
%         %     end
%         %     data = mean(HF,4);
%         data = permute(HF,[2 4 1 3]);
% 
%     else
%         data = reshape(data,fliplr(dsize));
% 
%         data = permute(data,[1,3,2]);
%         data = data(:,:,p);
%     end
% 
%     parm.dsize = dsize;
%     parm.ScanPt = ScanPt;
% end


% hObject    handle to pushbutton30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function con_dims_Callback(hObject, eventdata, handles)
% hObject    handle to con_dims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of con_dims as text
%        str2double(get(hObject,'String')) returns contents of con_dims as a double


% --- Executes during object creation, after setting all properties.
function con_dims_CreateFcn(hObject, eventdata, handles)
% hObject    handle to con_dims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton31.
function pushbutton31_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function filt_dim_Callback(hObject, eventdata, handles)
% hObject    handle to filt_dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filt_dim as text
%        str2double(get(hObject,'String')) returns contents of filt_dim as a double


% --- Executes during object creation, after setting all properties.
function filt_dim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filt_dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filt_cuts_Callback(hObject, eventdata, handles)
% hObject    handle to filt_cuts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filt_cuts as text
%        str2double(get(hObject,'String')) returns contents of filt_cuts as a double


% --- Executes during object creation, after setting all properties.
function filt_cuts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filt_cuts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in use_match.
function use_match_Callback(hObject, eventdata, handles)
% hObject    handle to use_match (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_match



function filt_match_Callback(hObject, eventdata, handles)
% hObject    handle to filt_match (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filt_match as text
%        str2double(get(hObject,'String')) returns contents of filt_match as a double


% --- Executes during object creation, after setting all properties.
function filt_match_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filt_match (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton33.
function pushbutton33_Callback(hObject, eventdata, handles)
P = handles.had_input.String;
P2 = evalin('base',P);
X = handles.had_output.String;
dims = size(P2);
switch length(dims)
    case 2
        H = hadamard(dims(2));
        n = dims(2);
        for i = 1:n
            X2(:,i) = mean(P2.*H(i,:),2);
        end
    case 3
        H = hadamard(dims(2));  
        n = dims(2);
        m = dims(3);
        for i = 1:n
            for j = 1:m
            X2(:,i,j) = mean(squeeze(P2(:,:,j)).*H(i,:),2);
            end
        end
    case 4
           H = hadamard(dims(3));
           n = dims(3);
           for i = 1:n
               K = H(i,:);
               K = permute(K,[1 3 2]);
               for j = 1:dims(4)

                   X2(:,:,i,j) = mean(P2(:,:,:,j).*K,3);
               end
               multiWaitbar('Converting',i/n);
           end
           multiWaitbar('CLOSEALL');
end
assignin('base',X,X2);


% hObject    handle to pushbutton33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function had_output_Callback(hObject, eventdata, handles)
% hObject    handle to had_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of had_output as text
%        str2double(get(hObject,'String')) returns contents of had_output as a double


% --- Executes during object creation, after setting all properties.
function had_output_CreateFcn(hObject, eventdata, handles)
% hObject    handle to had_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in hadamard.
function hadamard_Callback(hObject, eventdata, handles)
% hObject    handle to hadamard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~hObject.Value
    set(handles.text76,'String','dimensions');
else
    set(handles.text76,'String','group width');
    set(handles.planes,'Value',0);
end


% Hint: get(hObject,'Value') returns toggle state of hadamard


% --- Executes on button press in filter.
function filter_Callback(hObject, eventdata, handles)
% hObject    handle to filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filter


% --- Executes on button press in bband.
function bband_Callback(hObject, eventdata, handles)
set(handles.op_env,'Value',0);
% hObject    handle to bband (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bband



function bb_temp_Callback(hObject, eventdata, handles)
% hObject    handle to bb_temp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bb_temp as text
%        str2double(get(hObject,'String')) returns contents of bb_temp as a double


% --- Executes during object creation, after setting all properties.
function bb_temp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bb_temp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dbc.
function dbc_Callback(hObject, eventdata, handles)
% hObject    handle to dbc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dbc



function fs_Callback(hObject, eventdata, handles)
% hObject    handle to fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fs as text
%        str2double(get(hObject,'String')) returns contents of fs as a double


% --- Executes during object creation, after setting all properties.
function fs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in planes.
function planes_Callback(hObject, eventdata, handles)
set(handles.hadamard,'Value',0);
if ~hObject.Value
    set(handles.text76,'String','dimensions');
else
    set(handles.text76,'String','Theta  dTheta');
end
 
% hObject    handle to planes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of planes



function dtheta_Callback(hObject, eventdata, handles)
% hObject    handle to dtheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dtheta as text
%        str2double(get(hObject,'String')) returns contents of dtheta as a double


% --- Executes during object creation, after setting all properties.
function dtheta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dtheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ss_Callback(hObject, eventdata, handles)
% hObject    handle to ss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ss as text
%        str2double(get(hObject,'String')) returns contents of ss as a double


% --- Executes during object creation, after setting all properties.
function ss_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tp_m_Callback(hObject, eventdata, handles)
% hObject    handle to tp_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tp_m as text
%        str2double(get(hObject,'String')) returns contents of tp_m as a double


% --- Executes during object creation, after setting all properties.
function tp_m_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tp_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton35.
function pushbutton35_Callback(hObject, eventdata, handles)
X = evalin('base','X_c');
X = permute(X,[3 1 4 2]);
t = str2double(handles.tp_m.String);
if ndims(X) == 3
    M = X(:,:,t);
end
assignin('base','M',M);
% hObject    handle to pushbutton35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Venice.
function Venice_Callback(hObject, eventdata, handles)
% hObject    handle to Venice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
p = evalin('base','V');
pmean = mean(p(:));
p = p-pmean;
% if handles.op_env.Value
%     p = envelope(p')';
% end
% if handles.op_invert.Value
%     p = p*-1;
% end
if handles.dbc.Value
    p = real(20*log10(p./max(p(:))));
end
if handles.op_env.Value
    p = envelope(p')';
end
locx = [92 92 92 92 92 168 168 168 168 168 243 243 243 243 243];
locy = [50 100 150 200 250 50 100 150 200 250 50 100 150 200 250];
for i = 1:15
    x1 = locx(i);
    y1 = locy(i);
   % ymax1 = max(abs(p(:,x1)));
   ymax1 = p(y1,x1);
   ymax2 = ymax1;
   % ymax2 = max(abs(p(y1,:)));
    ymin2 = min(abs(p(y1,:)));
    ymin1 = min(abs(p(:,x1)));
    yhalf1 = (ymax1-ymin1)/2+ymin1;
    yhalf2 = (ymax2-ymin2)/2+ymin2;
    col = x1;
    row = y1;
    % axes(handles.axes1)
    % imagesc(p);
    pvert = p(:,col);
    phorz = p(row,:);
    pvert = resample(pvert,10,1);
    if handles.ext_fig.Value
        figure(26);
    else
        axes(handles.axes3)
    end
    plot(pvert,'b','LineWidth',2,'ButtonDownFcn',{@Plot_FFT,handles});
    if handles.op_fwhm_line.Value
        hold on;
        plot(ones(1,length(pvert))*yhalf1,'r--')
        hold off
    end
    phorz = resample(phorz,10,1);
    if handles.ext_fig.Value
        figure(27)
    else
        axes(handles.axes4)
    end
    hold off
    cla reset
    plot(phorz,'g','LineWidth',2,'ButtonDownFcn',{@Plot_FFT,handles});
    if handles.op_fwhm_line.Value
        hold on
        plot(ones(1,length(phorz))*yhalf2,'r--')
        hold off
    end
    %s2 = snum*10;
    v_first = pvert(1:row*10);
    v_second = pvert(row*10+1:end);
    h_first = phorz(1:col*10);
    h_second = phorz(col*10+1:end);
    x1 = find(v_first <= yhalf1,1,'last');
    x2 = find(v_second <= yhalf1,1)+row*10;
    % kgrid = evalin('base','kgrid');
    dy = str2double(handles.dy.String)/10;
    dx = str2double(handles.dx.String)/10;
    dz = str2double(handles.dz.String)/10;
    dt = str2double(handles.dt.String)/10;
    y1 = find(h_first <= yhalf2,1,'last');
    y2 = find(h_second <= yhalf2,1)+col*10;
    fwhm(i) = (x2-x1)*dy;
    fwhm2(i) = (y2-y1)*dx;
    maxval(i) = max([ymax1 ymax2]);
    
end
mets = [fwhm;fwhm2;maxval];
 assignin('base','peak',maxval);
 assignin('base','fwhm',fwhm);
 assignin('base','fwhm2',fwhm2);
assignin('base','mets',mets);



function fignum_Callback(hObject, eventdata, handles)
% hObject    handle to fignum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fignum as text
%        str2double(get(hObject,'String')) returns contents of fignum as a double


% --- Executes during object creation, after setting all properties.
function fignum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fignum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function figtits_Callback(hObject, eventdata, handles)
% hObject    handle to figtits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of figtits as text
%        str2double(get(hObject,'String')) returns contents of figtits as a double


% --- Executes during object creation, after setting all properties.
function figtits_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figtits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xlabel_Callback(hObject, eventdata, handles)
% hObject    handle to xlabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xlabel as text
%        str2double(get(hObject,'String')) returns contents of xlabel as a double


% --- Executes during object creation, after setting all properties.
function xlabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xlabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ylabel_Callback(hObject, eventdata, handles)
% hObject    handle to ylabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ylabel as text
%        str2double(get(hObject,'String')) returns contents of ylabel as a double


% --- Executes during object creation, after setting all properties.
function ylabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ylabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in keep_axes.
function keep_axes_Callback(hObject, eventdata, handles)
% hObject    handle to keep_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of keep_axes


% --- Executes on button press in add_cbar.
function add_cbar_Callback(hObject, eventdata, handles)
% hObject    handle to add_cbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of add_cbar


% --- Executes on button press in pushbutton37.
function pushbutton37_Callback(hObject, eventdata, handles)
fname = handles.savefigname.String;
if ~isempty(handles.fignum.String)
    num = str2double(handles.fignum.String);
    cb = handles.add_cbar.Value;
    ca = handles.convert.Value;
    if handles.keep_axes.Value == 0
        figure(num);
        a = gca;
        a.Visible = 'off';
    else
        figure(num)
        tits = handles.figtits.String;
        xlab  = handles.xlabel.String;
        ylab = handles.ylabel.String; 
        xlabel(xlab)
        ylabel(ylab)
        title(tits)
    end
    if cb
        colorbar;
    end
    %     if fname(end-2:end) == 'png'
    %         fname = [fname ' -transparent'];
    %         figure(num)
    %         set(gca,'Color','none');
    %     else
    %     end
    if ~isempty(handles.savefolder.String)
        pname = handles.savefolder.String;
        sname = [pname '\' fname];
    else
        sname = fname;
    end
    
    figure(num)
    set(gca,'Color','none')
    %  eval([ 'export_fig ' sname])
    export_fig(sname,'-transparent');
else
    disp('ERROR: You must select the figure number to export');
%     if ~isempty(handles.savefolder.String)
%         pname = handles.savefolder.String;
%         sname = [pname '\' fname];
%     else
%         sname = fname;
%     end
%     export_fig(sname) 
end
% hObject    handle to pushbutton37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in convert.
function convert_Callback(hObject, eventdata, handles)
% hObject    handle to convert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of convert


% --- Executes on button press in imshowbox.
function imshowbox_Callback(hObject, eventdata, handles)
% hObject    handle to imshowbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of imshowbox


% --- Executes on selection change in con_menu.
function con_menu_Callback(hObject, eventdata, handles)
v = handles.con_menu.String{handles.con_menu.Value};
switch v
    case 'Focus'
        set(handles.text76,'String','Focus dx')
          set(handles.had_input,'BackgroundColor',[1 1 1])
            set(handles.had_output,'BackgroundColor',[1 1 1])
    case 'Planes'
            set(handles.text76,'String','Theta dTheta')
            set(handles.had_input,'BackgroundColor',[.5 1 .5])
              set(handles.had_output,'BackgroundColor',[1 1 1])
    case 'Hadamard'
            set(handles.text76,'String','Element Width')
            set(handles.had_output,'BackgroundColor',[.5 1 .5])
            set(handles.had_input,'BackgroundColor',[1 1 1])
end
% hObject    handle to con_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns con_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from con_menu


% --- Executes during object creation, after setting all properties.
function con_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to con_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton38.
function pushbutton38_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
p = evalin('base','V');
pmean = mean(p(:));
p = p-pmean;
% if handles.op_env.Value
%     p = envelope(p')';
% end
% if handles.op_invert.Value
%     p = p*-1;
% end
if handles.dbc.Value
    p = real(20*log10(p./max(p(:))));
end
if handles.op_env.Value
    p = envelope(p')';
end
locx = [53 53 73 92   91  113 132 133 152 152 153 162];
locy = [100 40 60 140 100 40  141 131 100  85  70 171];
for i = 1:length(locx)
    x1 = locx(i);
    y1 = locy(i);
   % ymax1 = max(abs(p(:,x1)));
   ymax1 = p(y1,x1);
   ymax2 = ymax1;
   % ymax2 = max(abs(p(y1,:)));
    ymin2 = min(abs(p(y1,:)));
    ymin1 = min(abs(p(:,x1)));
    yhalf1 = (ymax1-ymin1)/2+ymin1;
    yhalf2 = (ymax2-ymin2)/2+ymin2;
    col = x1;
    row = y1;
    % axes(handles.axes1)
    % imagesc(p);
    pvert = p(:,col);
    phorz = p(row,:);
    pvert = resample(pvert,10,1);
    if handles.ext_fig.Value
        figure(26);
    else
        axes(handles.axes3)
    end
    plot(pvert,'b','LineWidth',2,'ButtonDownFcn',{@Plot_FFT,handles});
    if handles.op_fwhm_line.Value
        hold on;
        plot(ones(1,length(pvert))*yhalf1,'r--')
        hold off
    end
    phorz = resample(phorz,10,1);
    if handles.ext_fig.Value
        figure(27)
    else
        axes(handles.axes4)
    end
    hold off
    cla reset
    plot(phorz,'g','LineWidth',2,'ButtonDownFcn',{@Plot_FFT,handles});
    if handles.op_fwhm_line.Value
        hold on
        plot(ones(1,length(phorz))*yhalf2,'r--')
        hold off
    end
    %s2 = snum*10;
    v_first = pvert(1:row*10);
    v_second = pvert(row*10+1:end);
    h_first = phorz(1:col*10);
    h_second = phorz(col*10+1:end);
    x1 = find(v_first <= yhalf1,1,'last');
    x2 = find(v_second <= yhalf1,1)+row*10;
    % kgrid = evalin('base','kgrid');
    dy = str2double(handles.dy.String)/10;
    dx = str2double(handles.dx.String)/10;
    dz = str2double(handles.dz.String)/10;
    dt = str2double(handles.dt.String)/10;
    y1 = find(h_first <= yhalf2,1,'last');
    y2 = find(h_second <= yhalf2,1)+col*10;
    fwhm(i) = (x2-x1)*dy;
    fwhm2(i) = (y2-y1)*dx;
    maxval(i) = max([ymax1 ymax2]);
    
end
mets = [fwhm;fwhm2;maxval];
 assignin('base','peak',maxval);
 assignin('base','fwhm',fwhm);
 assignin('base','fwhm2',fwhm2);
assignin('base','mets',mets);



function had_input_Callback(hObject, eventdata, handles)
% hObject    handle to had_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of had_input as text
%        str2double(get(hObject,'String')) returns contents of had_input as a double


% --- Executes during object creation, after setting all properties.
function had_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to had_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in apply_mods.
function apply_mods_Callback(hObject, eventdata, handles)
% hObject    handle to apply_mods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of apply_mods


% --- Executes on button press in tr_delays.
function tr_delays_Callback(hObject, eventdata, handles)
% hObject    handle to tr_delays (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tr_delays


% --- Executes on button press in pushbutton39.
function pushbutton39_Callback(hObject, eventdata, handles)
M = evalin('base',handles.had_output.String); %bring in V(x,y,e)
kgrid = evalin('base','kgrid');
medium = evalin('base','medium');
source = evalin('base','source');
focus = str2num(handles.tr_focus.String); %bring in focus location
for i = 1:size(M,8) %envelope V in prep for delay determination
    M(:,i) = envelope(M(:,i));
end

[~, t] = max(M);
t2 = t-min(t);
%assignin('base','delays',t);

%Get Focal Delay Profile%

%                 for i = 1:length(element_index)
%                     steering_angle = atand((fpoint(2)-y0)/(element_index(i) + x0-fpoint(1)));
%                     burst1 = t0 + element_space * sin(steering_angle * pi/180)/...
%                         (ss * kgrid.dt);
%                     burst_offset(i) = burst1(i);
%                 end
%e_ind = element_index*(str2num(get(handles.trans_x_kerf,'String'))+width);

[x0 y0] = ind2sub(size(source.p_mask),find(source.p_mask));
xmid = median(x0);
x = x0-xmid;

w = focus(2)-y0; %depth in mm/10
h = (x-(focus(1)-xmid))*-1; %Height in mm/10
hyp = sqrt((h.*kgrid.dx).^2+(w.*kgrid.dy).^2); %hypoteneuse in mm/10
if length(medium.sound_speed) > 1
    medium.sound_speed = medium.sound_speed(1);
end
S = hyp./(medium.sound_speed(1)*kgrid.dt/1e4); %hypoteneuse in samples

burst_offset = hyp/medium.sound_speed(1)/kgrid.dt; %S;
burst_offset = abs(burst_offset-min(burst_offset))+1;
e = size(M,2);
if mod(e,2) == 1
    e = e+1;
end
b = burst_offset(round(linspace(1+e/2,length(h)-e/2,e)))';
delays = b-t2; %either this or b-t2
assignin('base','delays',delays);


% hObject    handle to pushbutton39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function tr_focus_Callback(hObject, eventdata, handles)
% hObject    handle to tr_focus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tr_focus as text
%        str2double(get(hObject,'String')) returns contents of tr_focus as a double


% --- Executes during object creation, after setting all properties.
function tr_focus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tr_focus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
