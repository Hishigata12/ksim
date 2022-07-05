function varargout = load_ct(varargin)
% LOAD_CT MATLAB code for load_ct.fig
%      LOAD_CT, by itself, creates a new LOAD_CT or raises the existing
%      singleton*.
%
%      H = LOAD_CT returns the handle to a new LOAD_CT or the handle to
%      the existing singleton*.
%
%      LOAD_CT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOAD_CT.M with the given input arguments.
%
%      LOAD_CT('Property','Value',...) creates a new LOAD_CT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before load_ct_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to load_ct_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help load_ct

% Last Modified by GUIDE v2.5 08-May-2022 23:46:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @load_ct_OpeningFcn, ...
                   'gui_OutputFcn',  @load_ct_OutputFcn, ...
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


% --- Executes just before load_ct is made visible.
function load_ct_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to load_ct (see VARARGIN)

% Choose default command line output for load_ct
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes load_ct wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = load_ct_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in cropct.
function cropct_Callback(hObject, eventdata, handles)
CT = evalin('base','CT');
CT2 = imcrop(CT);
assignin('base','CT2',CT2);
set(handles.cropped_box,'Value',1);
set(handles.pixelcount,'String',['Pixels: ' num2str(size(CT2))])
set(handles.xdims,'String',['X: ' num2str(size(CT2,2)/str2double(handles.x_range.String)) 'pix/mm']);
set(handles.ydims,'String',['Y: ' num2str(size(CT2,1)/str2double(handles.y_range.String)) 'pix/mm']);
% hObject    handle to cropct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in showorig.
function showorig_Callback(hObject, eventdata, handles)
CT = evalin('base','CT');
axes(handles.axes1);
imshow(CT);
% hObject    handle to showorig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in showcrop.
function showcrop_Callback(hObject, eventdata, handles)
CT = evalin('base','CT2');
axes(handles.axes1);
imshow(CT);
% figure;
% imshow(CT);
% hObject    handle to showcrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in loadct.
function loadct_Callback(hObject, eventdata, handles)
[file, path] = uigetfile('*.*');
CT = imread(fullfile(path,file));
if size(CT,3) > 1
    CT = rgb2gray(CT);
end
CT = im2double(CT);
assignin('base','CT',CT);
set(handles.pixelcount,'String',['Pixels: ' num2str(size(CT))])
set(handles.xdims,'String',['X: ' num2str(size(CT,2)/str2double(handles.x_range.String)) 'pix/mm']);
set(handles.ydims,'String',['Y: ' num2str(size(CT,1)/str2double(handles.y_range.String)) 'pix/mm']);
% hObject    handle to loadct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function x_range_Callback(hObject, eventdata, handles)
% hObject    handle to x_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_range as text
%        str2double(get(hObject,'String')) returns contents of x_range as a double


% --- Executes during object creation, after setting all properties.
function x_range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_range_Callback(hObject, eventdata, handles)
% hObject    handle to y_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_range as text
%        str2double(get(hObject,'String')) returns contents of y_range as a double


% --- Executes during object creation, after setting all properties.
function y_range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function topleft_Callback(hObject, eventdata, handles)
% hObject    handle to topleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of topleft as text
%        str2double(get(hObject,'String')) returns contents of topleft as a double


% --- Executes during object creation, after setting all properties.
function topleft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to topleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in set_range.
function set_range_Callback(hObject, eventdata, handles)
CT = evalin('base','CT');
set(handles.xdims,'String',['X: ' num2str(size(CT,2)/str2double(handles.x_range.String)) 'pix/mm']);
set(handles.ydims,'String',['Y: ' num2str(size(CT,1)/str2double(handles.y_range.String)) 'pix/mm']);
% hObject    handle to set_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
ss_min = str2num(get(handles.ss_min,'String'))*1000;
ss_max = str2num(get(handles.ss_max,'String'))*1000;
ss_t = str2num(get(handles.ss_threshold,'String'));
den_min = str2num(get(handles.density_min,'String'))*1000;
den_max = str2num(get(handles.density_max,'String'))*1000;
den_t = str2num(get(handles.density_threshold,'String'));
alpha_min = str2num(get(handles.alpha_min,'String'));
alpha_max = str2num(get(handles.alpha_max,'String'));
alpha_t = str2num(get(handles.alpha_threshold,'String'));
kgrid = evalin('base','kgrid');
ss_skull = str2num(get(handles.ss_skull,'String'))*1000;
den_skull = str2num(get(handles.den_skull,'String'))*1000;
alpha_skull = str2num(get(handles.alpha_skull,'String'));

%gets linear ranges for each value
ss_r = linspace(ss_t(1),ss_t(2),128);
alpha_r = linspace(alpha_t(1),alpha_t(2),128);
den_r = linspace(den_t(1),den_t(2),128);
ss = linspace(ss_min,ss_max,128);
alpha = linspace(alpha_min,alpha_max,128);
den = linspace(den_min,den_max,128);

if handles.use_3d.Value
    CT = evalin('base','ct3d');
    medium.sound_speed(1:kgrid.Nx,1:kgrid.Ny,1:kgrid.Nz) = ss_min;
    medium.alpha_coeff(1:kgrid.Nx,1:kgrid.Ny,1:kgrid.Nz) = alpha_min;
    medium.density(1:kgrid.Nx,1:kgrid.Ny,1:kgrid.Nz) = den_min;
    % medium.alpha_power(1:kgrid.Nx,1:kgrid.Ny) = medium.alpha_power(1,1);
    medium.alpha_power = 1.2;
    [Xq,Yq,Zq] = meshgrid(1:kgrid.Ny,1:kgrid.Nx,1:kgrid.Nz);
    CT2 = interp3(1:size(CT,2),1:size(CT,1),1:size(CT,3),double(CT),Xq,Yq,Zq);
    CT2 = CT2./max(CT2(:));
    for i = 1:size(CT2,1)
        for j = 1:size(CT2,2)
            for k = 1:size(CT2,3)
                if CT2(i,j,k) <= ss_t(1)
                    medium.sound_speed(i,j,k) = ss_min;
                elseif CT2(i,j,k) > ss_t(1) && CT2(i,j,k) < ss_t(2)
                    %             s = find(ss_r >= CT(i,j),1);
                    medium.sound_speed(i,j,k) = ss_max;
                else
                    medium.sound_speed(i,j,k) = ss_skull; %ss_max
                end
                if CT2(i,j,k) <= alpha_t(1)
                    medium.alpha_coeff(i,j,k) = alpha_min;
                    %         medium.alpha_power(t(1)+i,t(2)+j) = 1.05;
                elseif CT2(i,j,k) > alpha_t(1) && CT2(i,j,k) < alpha_t(2)
                    %             a = find(alpha_r >= CT(i,j),1);
                    medium.alpha_coeff(i,j,k) = alpha_max;
                    %             medium.alpha_power(t(1)+i,t(2)+j) = 1.10;
                else
                    medium.alpha_coeff(i,j,k) = alpha_skull;
                    %             medium.alpha_power(t(1)+i,t(2)+j) = 1.22;
                end
                if CT2(i,j,k) <= den_t(1)
                    medium.density(i,j,k) = den_min;
                elseif CT2(i,j,k) > den_t(1) && CT2(i,j,k) < den_t(2)
                    %             d = find(den_r >= CT(i,j),1);
                    medium.density(i,j,k) = den_max;
                else
                    medium.density(i,j,k) = den_skull; %den_max
                end
            end
        end
        multiWaitbar('Generating Medium',i/size(CT,1));
    end
else
    CT = evalin('base','CT2');
    t = str2num(handles.topleft.String);
    %Computes medium parameters
    % medium = evalin('base','medium');
    medium.sound_speed(1:kgrid.Nx,1:kgrid.Ny) = ss_min;
    medium.alpha_coeff(1:kgrid.Nx,1:kgrid.Ny) = alpha_min;
    medium.density(1:kgrid.Nx,1:kgrid.Ny) = den_min;
    % medium.alpha_power(1:kgrid.Nx,1:kgrid.Ny) = medium.alpha_power(1,1);
    medium.alpha_power = 1.2;

    %LINEAR SCALING REGIONS
    % for i = 1:size(CT,1)
    %     for j = 1:size(CT,2)
    %         if CT(i,j) <= ss_t(1)
    %             medium.sound_speed(t(1)+i,t(2)+ j) = ss_min;
    %         elseif CT(i,j) > ss_t(1) && CT(i,j) < ss_t(2)
    %             s = find(ss_r >= CT(i,j),1);
    %             medium.sound_speed(t(1)+i,t(2)+ j) = ss(s);
    %         else
    %             medium.sound_speed(t(1)+i,t(2)+ j) = 2440; %ss_max
    %         end
    %         if CT(i,j) <= alpha_t(1)
    %             medium.alpha_coeff(t(1)+i,t(2)+ j) = alpha_min;
    %         elseif CT(i,j) > alpha_t(1) && CT(i,j) < alpha_t(2)
    %             a = find(alpha_r >= CT(i,j),1);
    %             medium.alpha_coeff(t(1)+i,t(2)+ j) = alpha(a);
    %         else
    %             medium.alpha_coeff(t(1)+i,t(2)+ j) = alpha_max;
    %         end
    %         if CT(i,j) <= den_t(1)
    %             medium.density(t(1)+i,t(2)+ j) = den_min;
    %         elseif CT(i,j) > den_t(1) && CT(i,j) < den_t(2)
    %             d = find(den_r >= CT(i,j),1);
    %             medium.density(t(1)+i,t(2)+ j) = den(d);
    %         else
    %             medium.density(t(1)+i,t(2)+ j) = 1807; %den_max
    %         end
    %     end
    % end

    %STATIC MEAN REGIONS
    for i = 1:size(CT,1)
        for j = 1:size(CT,2)
            if CT(i,j) <= ss_t(1)
                medium.sound_speed(t(1)+i,t(2)+ j) = ss_min;
            elseif CT(i,j) > ss_t(1) && CT(i,j) < ss_t(2)
                %             s = find(ss_r >= CT(i,j),1);
                medium.sound_speed(t(1)+i,t(2)+ j) = ss_max;
            else
                medium.sound_speed(t(1)+i,t(2)+ j) = ss_skull; %ss_max
            end
            if CT(i,j) <= alpha_t(1)
                medium.alpha_coeff(t(1)+i,t(2)+ j) = alpha_min;
                %         medium.alpha_power(t(1)+i,t(2)+j) = 1.05;
            elseif CT(i,j) > alpha_t(1) && CT(i,j) < alpha_t(2)
                %             a = find(alpha_r >= CT(i,j),1);
                medium.alpha_coeff(t(1)+i,t(2)+ j) = alpha_max;
                %             medium.alpha_power(t(1)+i,t(2)+j) = 1.10;
            else
                medium.alpha_coeff(t(1)+i,t(2)+ j) = alpha_skull;
                %             medium.alpha_power(t(1)+i,t(2)+j) = 1.22;
            end
            if CT(i,j) <= den_t(1)
                medium.density(t(1)+i,t(2)+ j) = den_min;
            elseif CT(i,j) > den_t(1) && CT(i,j) < den_t(2)
                %             d = find(den_r >= CT(i,j),1);
                medium.density(t(1)+i,t(2)+ j) = den_max;
            else
                medium.density(t(1)+i,t(2)+ j) = den_skull; %den_max
            end
        end
    end
end
multiWaitbar('CLOSEALL')
assignin('base','medium',medium);
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function ss_min_Callback(hObject, eventdata, handles)
% hObject    handle to ss_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ss_min as text
%        str2double(get(hObject,'String')) returns contents of ss_min as a double


% --- Executes during object creation, after setting all properties.
function ss_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ss_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ss_max_Callback(hObject, eventdata, handles)
% hObject    handle to ss_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ss_max as text
%        str2double(get(hObject,'String')) returns contents of ss_max as a double


% --- Executes during object creation, after setting all properties.
function ss_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ss_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ss_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to ss_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ss_threshold as text
%        str2double(get(hObject,'String')) returns contents of ss_threshold as a double


% --- Executes during object creation, after setting all properties.
function ss_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ss_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function density_min_Callback(hObject, eventdata, handles)
% hObject    handle to density_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of density_min as text
%        str2double(get(hObject,'String')) returns contents of density_min as a double


% --- Executes during object creation, after setting all properties.
function density_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to density_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function density_max_Callback(hObject, eventdata, handles)
% hObject    handle to density_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of density_max as text
%        str2double(get(hObject,'String')) returns contents of density_max as a double


% --- Executes during object creation, after setting all properties.
function density_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to density_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function density_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to density_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of density_threshold as text
%        str2double(get(hObject,'String')) returns contents of density_threshold as a double


% --- Executes during object creation, after setting all properties.
function density_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to density_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alpha_min_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha_min as text
%        str2double(get(hObject,'String')) returns contents of alpha_min as a double


% --- Executes during object creation, after setting all properties.
function alpha_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alpha_max_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha_max as text
%        str2double(get(hObject,'String')) returns contents of alpha_max as a double


% --- Executes during object creation, after setting all properties.
function alpha_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alpha_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha_threshold as text
%        str2double(get(hObject,'String')) returns contents of alpha_threshold as a double


% --- Executes during object creation, after setting all properties.
function alpha_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cropped_box.
function cropped_box_Callback(hObject, eventdata, handles)
if hObject.Value
    CT = evalin('base','CT2');
else 
    CT = evalin('base','CT');
end
set(handles.pixelcount,'String',['Pixels: ' num2str(size(CT))])
set(handles.xdims,'String',['X: ' num2str(size(CT,2)/str2double(handles.x_range.String)) 'pix/mm']);
set(handles.ydims,'String',['Y: ' num2str(size(CT,1)/str2double(handles.y_range.String)) 'pix/mm']);

% hObject    handle to cropped_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cropped_box



function intx_Callback(hObject, eventdata, handles)
% hObject    handle to intx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of intx as text
%        str2double(get(hObject,'String')) returns contents of intx as a double


% --- Executes during object creation, after setting all properties.
function intx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to intx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function inty_Callback(hObject, eventdata, handles)
% hObject    handle to inty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inty as text
%        str2double(get(hObject,'String')) returns contents of inty as a double


% --- Executes during object creation, after setting all properties.
function inty_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in int_button.
function int_button_Callback(hObject, eventdata, handles)
if handles.cropped_box.Value
    CT = evalin('base','CT2');
else 
    CT = evalin('base','CT');
end

in = [str2double(handles.inty.String) str2double(handles.intx.String)];
dims = size(CT);
[x, z] = meshgrid(1:dims(1),1:dims(2));
inval = dims./in;
% [x2,z2] = meshgrid(1:(1/in(1)):dims(1),1:(1/in(2)):dims(2));
[x2,z2] = meshgrid(1:inval(1):dims(1),1:inval(2):dims(2));
CT2 = interp2(x,z,CT',x2,z2);
CT2 = rot90(CT2,3);
CT2 = fliplr(CT2);
         
if handles.cropped_box.Value
    assignin('base','CT2',CT2)
else 
    assignin('base','CT',CT2);
end
set(handles.pixelcount,'String',['Pixels: ' num2str(size(CT2))])
set(handles.xdims,'String',['X: ' num2str(size(CT2,2)/str2double(handles.x_range.String)) 'pix/mm']);
set(handles.ydims,'String',['Y: ' num2str(size(CT2,1)/str2double(handles.y_range.String)) 'pix/mm']);
% hObject    handle to int_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in fliplr_button.
function fliplr_button_Callback(hObject, eventdata, handles)
if handles.cropped_box.Value
    CT = evalin('base','CT2');
else 
    CT = evalin('base','CT');
end

CT = fliplr(CT);

if handles.cropped_box.Value
    assignin('base','CT2',CT)
else 
    assignin('base','CT',CT);
end
set(handles.pixelcount,'String',['Pixels: ' num2str(size(CT))])
set(handles.xdims,'String',['X: ' num2str(size(CT,2)/str2double(handles.x_range.String)) 'pix/mm']);
set(handles.ydims,'String',['Y: ' num2str(size(CT,1)/str2double(handles.y_range.String)) 'pix/mm']);

% hObject    handle to fliplr_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in flipud_button.
function flipud_button_Callback(hObject, eventdata, handles)
if handles.cropped_box.Value
    CT = evalin('base','CT2');
else 
    CT = evalin('base','CT');
end

CT = flipud(CT);

if handles.cropped_box.Value
    assignin('base','CT2',CT)
else 
    assignin('base','CT',CT);
end
set(handles.pixelcount,'String',['Pixels: ' num2str(size(CT))])
set(handles.xdims,'String',['X: ' num2str(size(CT,2)/str2double(handles.x_range.String)) 'pix/mm']);
set(handles.ydims,'String',['Y: ' num2str(size(CT,1)/str2double(handles.y_range.String)) 'pix/mm']);
% hObject    handle to flipud_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function ss_skull_Callback(hObject, eventdata, handles)
% hObject    handle to ss_skull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ss_skull as text
%        str2double(get(hObject,'String')) returns contents of ss_skull as a double


% --- Executes during object creation, after setting all properties.
function ss_skull_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ss_skull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function den_skull_Callback(hObject, eventdata, handles)
% hObject    handle to den_skull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of den_skull as text
%        str2double(get(hObject,'String')) returns contents of den_skull as a double


% --- Executes during object creation, after setting all properties.
function den_skull_CreateFcn(hObject, eventdata, handles)
% hObject    handle to den_skull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alpha_skull_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_skull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha_skull as text
%        str2double(get(hObject,'String')) returns contents of alpha_skull as a double


% --- Executes during object creation, after setting all properties.
function alpha_skull_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_skull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rot_Callback(hObject, eventdata, handles)
% hObject    handle to rot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rot as text
%        str2double(get(hObject,'String')) returns contents of rot as a double


% --- Executes during object creation, after setting all properties.
function rot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rotbutt.
function rotbutt_Callback(hObject, eventdata, handles)
if handles.cropped_box.Value
    CT = evalin('base','CT2');
else
CT = evalin('base','CT');
end
rot = str2double(handles.rot.String);
CT2 = imrotate(CT,rot,'bilinear','crop');
assignin('base','CT2',CT2);

% hObject    handle to rotbutt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

        


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
d = handles.dicom_folder.String{1};
n = handles.dicom_fname.String{1};
f = fullfile(d,n);
cd(d);
a = dir;
for i = 1:length(a)-3
     filename = sprintf([n '_%03d.dcm'], i);
   X(:,:,i) = dicomread(filename);
end
   
assignin('base','ct3d',X)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function dicom_fname_Callback(hObject, eventdata, handles)
% hObject    handle to dicom_fname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dicom_fname as text
%        str2double(get(hObject,'String')) returns contents of dicom_fname as a double


% --- Executes during object creation, after setting all properties.
function dicom_fname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dicom_fname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dicom_folder_Callback(hObject, eventdata, handles)
% hObject    handle to dicom_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dicom_folder as text
%        str2double(get(hObject,'String')) returns contents of dicom_folder as a double


% --- Executes during object creation, after setting all properties.
function dicom_folder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dicom_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in use_3d.
function use_3d_Callback(hObject, eventdata, handles)
% hObject    handle to use_3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_3d
