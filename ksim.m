
function varargout = ksim(varargin)
% KSIM MATLAB code for ksim.fig
%      KSIM, by itself, creates a new KSIM or raises the existing
%      singleton*.
%
%      H = KSIM returns the handle to a new KSIM or the handle to
%      the existing singleton*.
%
%      KSIM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KSIM.M with the given input arguments.
%
%      KSIM('Property','Value',...) creates a new KSIM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ksim_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ksim_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ksim

% Last Modified by GUIDE v2.5 07-Jun-2022 16:59:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ksim_OpeningFcn, ...
                   'gui_OutputFcn',  @ksim_OutputFcn, ...
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

% TR = [22.84 22.6 22.4 22.2 22.04 21.88 21.76 21.64 21.56 21.48 21.44 21.4 21.36 21.4 21.44 21.48 21.56 21.64 21.76 21.88 22.04 22.2 22.4 22.6 22.84];
% --- Executes just before ksim is made visible.
function ksim_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ksim (see VARARGIN)

% Choose default command line output for ksim
handles.output = hObject;
if ismember('sensor_data',evalin('base','who'))
    sensor_data = evalin('base','sensor_data');
    set(handles.data_menu,'String',fieldnames(sensor_data));
    if isstruct(sensor_data)
        sensor_data = sensor_data.p;
    end
    update_sensormenu(hObject, eventdata, handles, sensor_data);
    if ismember('kgrid',evalin('base','who'))
        kgrid = evalin('base','kgrid');
        update_text(hObject, eventdata, handles, sensor_data, kgrid);
    end
end

set(handles.medium_layers,'Visible','off')
set(handles.trans_focus,'Visible','on');
set(handles.trans_focus,'String','200 270')
set(handles.sensormenu,'String',1:400*400);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ksim wait for user response (see UIRESUME)
% uiwait(handles.figure1);
function update_text(hObject, eventdata, handles, sensor_data, kgrid);
sensor_size = size(sensor_data);
% new_size = [num2str(sensor_size(1)), ' Sensors ', num2str(sensor_size(2)),' Points'];
% set(handles.time_sensor_text,'String',new_size);

dt = kgrid.dt*1e6; %dt in [us]
set(handles.dt_text,'String',[num2str(dt), ' usec/pt']);

function update_sensormenu(hObject, eventdata, handles, sensor_data);
S = linspace(1,size(sensor_data,1),size(sensor_data,1)); %Sensor numbers
set(handles.sensormenu,'String',S);


% --- Outputs from this function are returned to the command line.
function varargout = ksim_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in sensormenu.
function sensormenu_Callback(hObject, eventdata, handles)
% hObject    handle to sensormenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sensormenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sensormenu


% --- Executes during object creation, after setting all properties.
function sensormenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensormenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in show_sensor.
function show_sensor_Callback(hObject, eventdata, handles)
kgrid = evalin('base','kgrid');
%sensor_data = evalin('base','sensor_data');
val = handles.data_menu.String{handles.data_menu.Value};

%Gets vars depending on parameter
switch val
    case 'p'
        pressure = evalin('base','pressure');
     %   clims = [max(abs(sensor_data(:)))*2*-1, max(abs(sensor_data(:)))*2];%[-1 1];
        ylab = 'Pressure';
        ylab1 = 'Sensor Number';
        xlab1 = 'Time Step';
        tits1 = 'Pressure';
        plot1 = 0;
    case 'ux'
        sensor_data = sensor_data.ux;
      %  clims = [max(abs(sensor_data(:)))*2*-1, max(abs(sensor_data(:)))*2];
         ylab = 'Velocity X';
          ylab1 = 'Sensor Number';
        xlab1 = 'Time Step';
        tits1 = 'X Velocity';
         plot1 = 0;
    case 'uy'
        sensor_data = sensor_data.uy;
      %  clims = [max(abs(sensor_data(:)))*2*-1, max(abs(sensor_data(:)))*2];
         ylab = 'Velocity Y';
           ylab1 = 'Sensor Number';
        xlab1 = 'Time Step';
        tits1 = 'Y Velocity';
         plot1 = 0;
    case 'p_final'
        sensor_data = sensor_data.p_final;
      %  clims = [max(abs(sensor_data(:)))*2*-1, max(abs(sensor_data(:)))*2];
        ylab1 = 'X Position';
        xlab1 = 'Y Position';
        tits1 = 'Final Pressure';
         plot1 = 1;
    case 'ux_final'
        sensor_data = sensor_data.ux_final;
       % clims = [max(abs(sensor_data(:)))*2*-1, max(abs(sensor_data(:)))*2];
        ylab1 = 'X Position';
        xlab1 = 'Y Position';
        tits1 = 'Final Velocity X';
         plot1 = 1;
    case 'uy_final'
        sensor_data = sensor_data.uy_final;
      %  clims = [max(abs(sensor_data(:)))*2*-1, max(abs(sensor_data(:)))*2];
        ylab1 = 'X Position';
        xlab1 = 'Y Position';
        tits1 = 'Final Velocity Y';
        plot1 = 1;
    case 'p_max'
        sensor_data = sensor_data.p_max;
        xlab1 = 'Sensor Number';
        ylab1 = '';
        tits1 = 'Max Pressure';
        plot1 = 2;
    case 'p_max_all'
        sensor_data = sensor_data.uy_final;
        tits1 = 'Max Pressure';
        ylab1 = 'X Position';
        xlab1 = 'Y Position';
        plot1 = 1;
    case 'Ix'
            sensor_data = sensor_data.Ix;
     %   clims = [max(abs(sensor_data(:)))*2*-1, max(abs(sensor_data(:)))*2];%[-1 1];
        ylab = 'Intensity';
        ylab1 = 'Sensor Number';
        xlab1 = 'Time Step';
        tits1 = 'Intensity X';
        plot1 = 0;
    case 'Iy'
            sensor_data = sensor_data.Iy;
     %   clims = [max(abs(sensor_data(:)))*2*-1, max(abs(sensor_data(:)))*2];%[-1 1];
        ylab = 'Intensity';
        ylab1 = 'Sensor Number';
        xlab1 = 'Time Step';
        tits1 = 'Intensity Y';
        plot1 = 0;
end
sensor = evalin('base','sensor');
% if strcmp('Custom',sensor.type)
    sensor.mask = flipud(sensor.mask);
% end
    

% Gets CLims
if isempty(get(handles.display_clims,'String'))
    clims = [max(abs(sensor_data(:)))*2*-1, max(abs(sensor_data(:)))*2];
else
    clims = str2num(handles.display_clims.String);
end

%Plots on Axes 1
% axes(handles.axes1)
% if plot1 == 0
% imagesc(kgrid.t_array,1:size(sensor_data,1),sensor_data,clims);
% elseif plot1 == 1
%     imagesc(sensor_data,clims);
% elseif plot1 == 2
%     if handles.twodsensor.Value
%         sensor_data = reshape(sensor_data,[sqrt(length(sensor_data)), sqrt(length(sensor_data))]);
%         imagesc(sensor_data);
%     else
%     plot(sensor_data)
%     end
% end
% colormap('jet');
% ylabel(ylab1);    
% xlabel(xlab1);
% title(tits1);
% colorbar;

%Plots on axes 2
if plot1 == 0
axes(handles.axes2)
Sn = round(get(handles.sensormenu,'Value'));
lin_pres = reshape(pressure,[size(pressure,1)*size(pressure,2) size(pressure,3)]);
plot(kgrid.t_array,lin_pres(Sn,:))
%plot(kgrid.t_array,sensor_data(Sn,:))
ylabel(ylab);
xlabel('Time');
title('Pressure over time')
end

%Plots on axes 3
invertp3 = 1;

if invertp3
    X1 = 2;
    Y1 = 1;
else
    X1 = 1;
    Y1 = 2;
end
if plot1 == 0
sensor = evalin('base','sensor');
sensor.mask = sensor.mask*100;
if handles.ef.Value
    figure(6);
else
    axes(handles.axes3);
end
hold off
sensor.mask(2,:) = sensor.mask(2,:) - min(sensor.mask(2,:));
scatter(sensor.mask(X1,:),sensor.mask(Y1,:),'k');
title('Sensor Location')
hold all;
% if strcmp(sensor.type,'Circle') || strcmp(sensor.type,'Cart_Circle')
%     shift = find(sensor.mask(1,:) == min(sensor.mask(1,:)));
%     sensor.mask = circshift(sensor.mask,-shift+1,2);
% end
scatter(sensor.mask(X1, Sn),sensor.mask(Y1,Sn),'r','filled');
% scatter(sensor.mask(1, Sn + shift - 1),sensor.mask(2,Sn + shift - 1),'r','filled');
xsize = [min(kgrid.x_vec)*100,max(kgrid.x_vec)*100];
ysize = [min(kgrid.y_vec)*100,max(kgrid.y_vec)*100];
ysize = ysize - min(ysize);
xlim(ysize); ylim(xsize);
colorbar;
end


% hObject    handle to show_sensor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in sensor.
function sensor_Callback(hObject, eventdata, handles)
sensor = evalin('base','sensor');
kgrid = evalin('base','kgrid');
if kgrid.Nz == 0
    axes(handles.axes2);
    scatter(sensor.mask(2,:)-min(sensor.mask(2,:)),sensor.mask(1,:));
    ysize = [min(kgrid.x_vec),max(kgrid.x_vec)];
    xsize = [min(kgrid.y_vec),max(kgrid.y_vec)];
    xlim(xsize-(min(xsize))); ylim(ysize);
else
    sizes = size(kgrid.x);
    m = zeros(sizes);
    if numel(sensor.mask) ~= prod(sizes)
        %find the sensor plane
        f(1) = numel(unique(sensor.mask(1,:)));
        f(2) = numel(unique(sensor.mask(2,:)));
        f(3) = numel(unique(sensor.mask(3,:)));
        if f(1) == 1
              p = find(kgrid.x_vec == sensor.mask(1,1));
            m(p,:,:) = ones(sizes(2),sizes(3));
        elseif f(2) == 1
              p = find(kgrid.y_vec == sensor.mask(2,1));
            m(:,p,:) = ones(sizes(1),sizes(3));
        elseif f(3) == 1
            p = find(kgrid.z_vec == sensor.mask(3,1));
            m(:,:,p) = ones(sizes(1),sizes(2));
            %             for i = 1:sizes(1)

        end
           voxelPlot(m);
    end
end


            % hObject    handle to sensor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in source.
function source_Callback(hObject, eventdata, handles)
source = evalin('base','source');
kgrid = evalin('base','kgrid');
if ~ismember('p0',fieldnames(source))
    if kgrid.Nz == 0
        source.p0 = zeros(kgrid.Nx,kgrid.Ny);
    else 
        source.p0 = zeros(kgrid.Nx,kgrid.Ny,kgrid.Nz);
    end
end
if ismember('p_mask',fieldnames(source))
    source.p0 = source.p0+source.p_mask;
end
xsize = ([min(kgrid.x_vec),max(kgrid.x_vec)]-min(kgrid.x_vec))*100;
ysize = ([min(kgrid.y_vec),max(kgrid.y_vec)]-min(kgrid.y_vec))*100;
x = (kgrid.x_vec-min(kgrid.x_vec))*100;
y = (kgrid.y_vec-min(kgrid.y_vec))*100;
if ndims(source.p_mask) == 2
    if handles.ef.Value
        figure(2);
        X = ones(1,size(source.p,1)).*str2double(handles.trans_y0.String);
        Y = linspace(0,1,str2double(handles.trans_elements.String));
        Y = Y.*(str2double(handles.trans_x_kerf.String)+str2double(handles.trans_width.String)).*str2double(handles.trans_elements.String)+1;
        %         Y = 1:str2double(handles.trans_x_kerf.String):str2double(handles.trans_elements.String);
        Y = Y+str2double(handles.trans_x0.String)-median(Y);
        X = X.*kgrid.dy*100;
        Y = Y.*kgrid.dx*100;
        scatter(X,Y,'filled','s','MarkerFaceColor','b')
        % imagesc(y,x,source.p0)
        %         scatter(X,Y,'filled','MarkerEdgeColor',[0.8 0.2 0.2],...
        %               'MarkerFaceColor',[0.8 0.2 0.2],'square');
        xlim(round(ysize)); ylim(round(xsize));
    else
        axes(handles.axes2);
        imagesc(y,x,source.p0);
        xlim(ysize); ylim(xsize);
    end

    ylabel('Lateral cm')
    xlabel('Depth cm')
    colormap([1 1 1;1 0.2 0.2])
    title('Source')
else
    voxelPlot(source.p0);
end
if ismember('p',fieldnames(source))
    axes(handles.axes3)
    hold off;

    % trans = handles.trans_menu.String{handles.trans_menu.Value};

    axes(handles.axes3)
    yw = str2double(handles.trans_width.String);
    zw = str2double(handles.trans_length.String);
    shift = linspace(0,1,size(source.p,1)/yw/zw+1);
    if handles.ef.Value
        img = figure;

        for i = 1:size(source.p,1)/yw/zw
            pos = [0 1-shift(i+1) 1 1/size(source.p,1)*yw*zw];
            subplot('Position',pos)
            if length(kgrid.t_array) > size(source.p,2)
                img = plot(kgrid.t_array(1:size(source.p,2)),source.p(1+(i-1)*yw*zw,:),'k','LineWidth',3);
            else
                img = plot(kgrid.t_array,source.p(1+(i-1)*yw*zw,1:length(kgrid.t_array)),'k','LineWidth',3);
            end
            xticks([]);
            yticks([]);
            set(gca,'Color','none')
            xlabel(''); ylabel('');
        end
        % %
        %         if length(kgrid.t_array) > size(source.p,2)
        %             plot(kgrid.t_array(1:size(source.p,2))*1e6,source.p)
        %         else
        %             plot(kgrid.t_array*1e6,source.p)
        %         end

        xlabel('Time [us]')
        ylabel('Magnitude')
        title('Excitation')
    end



end
% hObject    handle to source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in simulate.
function simulate_Callback(hObject, eventdata, handles)
%import and convert variables from GUI
medium = evalin('base','medium');
sensor = evalin('base','sensor');
kgrid = evalin('base','kgrid');
rec = get(handles.record_movie,'Value'); %Changed from Boolean
mesh = get(handles.sim_mesh,'Value');
mask1 = get(handles.sim_mask,'Value');
% moviename = get(handles.rec_name,'String');
moviename = 'poopooman';
if mask1
    mask = 'on';
else
    mask = 'off';
end
pml = get(handles.sim_pml,'Value');
if ~isempty(handles.sim_scale.String)
    scale = str2num(get(handles.sim_scale,'String'));
else
    scale = 'auto';
end
freq = str2double(get(handles.sim_freq,'String'));
fps = str2double(get(handles.sim_fps,'String'));
psim = get(handles.sim_psim,'Value');

if handles.use_gpu.Value
    dc = 'gpuArray-single';
else 
    dc = 'off';
end

%Prep vars for use
fdir = get(handles.record_dir,'String');
fname = get(handles.record_filename,'String');
fname2 = fullfile(fdir,fname);
rec_args = {'RecordMovie',rec,'MovieName',fname2,'PlotScale',scale,'PlotFreq',freq,...
    'MovieProfile','MPEG-4','MovieArgs',{'FrameRate',fps}};
mesh_args = {'MeshPlot', mesh};
mask_args = {'DisplayMask',mask};
pml_args = {'PlotPML',pml};
plot_args = {'PlotSim',psim};
dc_args = {'DataCast',dc};
% input_args = {'PlotScale',scale,'PlotFreq',freq};
% input_args = [rec_args mesh_args, mask_args, pml_args];
% kgrid.dt = str2double(get(handles.dt_sim,'String'))*1e-6;
if handles.sim_3D.Value
    input_args = [rec_args, mask_args, pml_args, dc_args, plot_args];
    %     if ismember('transducer',evalin('base','who'))
    if handles.use_trans.Value
        transducer = evalin('base','transducer');
        [sensor_data] = kspaceFirstOrder3D(kgrid, medium, transducer, sensor, input_args{:});
    else
        source = evalin('base','source');
        sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
    end
else
    input_args = [rec_args mesh_args, mask_args, pml_args, dc_args, plot_args];
    source = evalin('base','source');
    sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
end
axes(handles.axes1);
% sensor_data.p2 = ones(size(sensor_data.p));
% if ~handles.sim_3D.Value
    if handles.use_gpu.Value
        sensor_data.p = double(gather(sensor_data.p));
    end
%     else
%         p2 = sensor_data.p;
%     end
% for i = 1:size(pressure,1)
%     for j = 1:size(pressure,2)
%         p_max(i,j) = max(abs(pressure(i,j,:)));
%     end
% end

%     for i = 1:size(p2,3)
%         sensor_data.p2(:,:,i) = envelope(p2(:,:,i));
%     end
% end
  %  assignin('base','sensor_data',sensor_data);
    set(handles.data_menu,'String',fieldnames(sensor_data));
if isstruct(sensor_data)
    sensor_data = sensor_data.p;
end

update_sensormenu(hObject, eventdata, handles, sensor_data);

%update_text(hObject, eventdata, handles, sensor_data, kgrid);
set(handles.tslide,'Max',size(sensor_data,2));
set(handles.tslide,'Min',1);
set(handles.tslide,'SliderStep',[1/size(sensor_data,2) 0.1]);
if strcmp(handles.wave_type.String{handles.wave_type.Value},'Focus')
    focus = str2num(handles.trans_focus.String);
    sRow = focus(1)/(kgrid.dx*1e4);
    sCol = round(focus(2)/(kgrid.dy*1e4));
    SNUM = kgrid.Nx/(kgrid.dx*1e4)*(sCol-1)+sRow;
    if handles.sim_3D.Value
       % if handles.all_sensor.Value
       if sensor.all
            SNUM2 = sub2ind([kgrid.Ny kgrid.Nx kgrid.Nz],focus(2), focus(1), focus(3));
       else
           plane = sensor.plane;
           switch plane
               case 'XY'
                    SNUM2 = sub2ind([kgrid.Ny kgrid.Nx],focus(2), focus(1));
               case 'XZ'
                     SNUM2 = sub2ind([kgrid.Nz kgrid.Nx],focus(3), focus(1));
               case 'YZ'
                    SNUM2 = sub2ind([kgrid.Ny kgrid.Nz],focus(2), focus(3));
           end
%            if handles.p_menu.Value == 1 %XY
%                if strcmp(sensor.plane,'XY')
%                    SNUM2 = sub2ind([kgrid.Ny kgrid.Nx],focus(2), focus(1));
%                elseif handles.p_menu.Value == 2 %XZ
%                    SNUM2 = sub2ind([kgrid.Nz kgrid.Nx],focus(3), focus(1));
%                elseif handles.p_menu.Value == 3 %YZ
%                    SNUM2 = sub2ind([kgrid.Ny kgrid.Nz],focus(2), focus(3));
%                end
        end
    else
    SNUM2 = sub2ind([kgrid.Nx kgrid.Ny],focus(1), focus(2));
    end
    set(handles.sensormenu,'Value',SNUM2)
end
if ~handles.sim_3D.Value
    pressure = reshape(sensor_data,[kgrid.Nx,kgrid.Ny,kgrid.Nt]);
else
   % if handles.all_sensor.Value
   plane = sensor.plane;
   if sensor.all
       pressure = reshape(sensor_data,[kgrid.Nx,kgrid.Ny,kgrid.Nz,kgrid.Nt]);
       pressure = permute(pressure,[2 1 3 4]);
   else
       switch plane
           case 'XY'
               pressure = reshape(sensor_data,[kgrid.Nx,kgrid.Ny,kgrid.Nt]);
               pressure = permute(pressure,[2 1 3]);
           case 'XZ'
               pressure = reshape(sensor_data,[kgrid.Nx,kgrid.Nz,kgrid.Nt]);
               pressure = permute(pressure,[2 1 3]);
           case 'YZ'
               pressure = reshape(sensor_data,[kgrid.Ny,kgrid.Nz,kgrid.Nt]);
       end

       %         if handles.p_menu.Value == 1 %XY
       %             pressure = reshape(sensor_data,[kgrid.Nx,kgrid.Ny,kgrid.Nt]);
       %             pressure = permute(pressure,[2 1 3]);
       %         elseif handles.p_menu.Value == 2 %XZ
       %             pressure = reshape(sensor_data,[kgrid.Nx,kgrid.Nz,kgrid.Nt]);
       %             pressure = permute(pressure,[2 1 3]);
       %         elseif handles.p_menu.Value == 3 %YZ
       %             pressure = reshape(sensor_data,[kgrid.Ny,kgrid.Nz,kgrid.Nt]);
       %         end
   end
end

%         if numel(sensor_data) == kgrid.Nx*kgrid.Ny*kgrid.Nt
%             pressure = reshape(sensor_data,[kgrid.Nx,kgrid.Ny,kgrid.Nt]);
%             clear sensor_data
%         end
%         if handles.sim_3D.Value
%             pressure = permute(pressure,[2 1 3 4]);
%         end
assignin('base','pressure',pressure);

% hObject    handle to simulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function tp_Callback(hObject, eventdata, handles)
% hObject    handle to tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tp as text
%        str2double(get(hObject,'String')) returns contents of tp as a double


% --- Executes during object creation, after setting all properties.
function tp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in show_tp.
function show_tp_Callback(hObject, eventdata, handles)
kgrid = evalin('base','kgrid');
sensor_data = evalin('base','sensor_data');
sensor = evalin('base','sensor');
if isstruct(sensor_data)
    if handles.show_envelope.Value
        sensor_data = imgaussfilt(abs(sensor_data.p),0.5);
    else
        sensor_data = sensor_data.p;
    end
end

axes(handles.axes1)
if handles.twodsensor.Value
    v = sqrt(size(sensor_data,1));
    sensor_data2 = reshape(sensor_data,[v,v,size(sensor_data,2)]);
    imagesc(sensor_data2(:,:,1));
else
    imagesc(kgrid.t_array.*1e6,1:size(sensor_data,1),sensor_data, [-1, 1]);
    colormap('hotcold');
    ylabel('Sensor Number');
    xlabel('Time Step');
    colorbar;
end


Tn = str2num(get(handles.tp,'String'));
ymax1 = 0;
ymin1 = 0;

if size(sensor.mask,2) > kgrid.Ny
    xline = ((str2double(handles.sensor_2.String)-1)*kgrid.Ny+1):(str2double(handles.sensor_2.String))*kgrid.Ny;
    yline = str2double(handles.sensor_1.String):kgrid.Nx:(str2double(handles.sensor_1.String)+(kgrid.Nx*kgrid.Ny-1));
else
    xline = 1:kgrid.Ny;
    yline = 0;
end

for i = 1:length(Tn)
    Tp = Tn(i);
    if handles.ef.Value
        figure(5);
    else
        axes(handles.axes2)
    end
    if i == 1
        if handles.yline.Value
            plot(1:length(yline),sensor_data(yline,Tp),'k','LineWidth',2.5)
        else
            plot(1:length(xline),sensor_data(xline,Tp),'k','LineWidth',2.5)
        end
    else
        if handles.yline.Value
             handles.axes2.Children.YData = sensor_data(yline,Tp);
        else
            handles.axes2.Children.YData = sensor_data(xline,Tp);
        end
    end
    t = num2str(kgrid.t_array(Tp).*1e6);
    title(['tp = ', t, ' us']);
    ylabel('Pressure');
    xlabel('Sensor');
    if handles.holdmax.Value
        if handles.max_global.Value
            ymin1 = min(sensor_data(:));
            ymax1 = max(sensor_data(:));
        else
            ymax = max(sensor_data(:,Tp));
            ymin = min(sensor_data(:,Tp));
            if ymax > ymax1
                ymax1 = ymax;
            end
            if ymin < ymin1
                ymin1 = ymin;
            end
        end
        ylim([ymin1 ymax1]);
    end
    
    sensor = evalin('base','sensor');
    if strcmp('Custom',sensor.type)
        sensor.mask = flipud(sensor.mask);
    end
    
    if handles.ef.Value
        figure(55);
    else
        axes(handles.axes3);
    end
    hold off
    color_vals = (sensor_data(:,Tp)-min(sensor_data(:)))./max(sensor_data(:));
    color_vals = sensor_data(:,Tp);
    c = [color_vals, zeros(length(color_vals),2)];
    if i ==1
        scatter(sensor.mask(1,:),sensor.mask(2,:),[],color_vals);
        colormap('hotcold');
    else
        handles.axes3.Children.CData = color_vals;
    end
    hold all;
    shift = find(sensor.mask(1,:) == min(sensor.mask(1,:)));
    % scatter(sensor.mask(1,shift - 1),sensor.mask(2, shift - 1),'r','filled');
    xsize = [min(kgrid.x_vec),max(kgrid.x_vec)];
    ysize = [min(kgrid.y_vec),max(kgrid.y_vec)];
    xlim(xsize); ylim(ysize);
    
    if handles.twodsensor.Value
    axes(handles.axes1);
    imagesc(sensor_data2(:,:,Tp)); 
    end
    
end

    

% hObject    handle to show_tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function dt_sim_Callback(hObject, eventdata, handles)
% hObject    handle to dt_sim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dt_sim as text
%        str2double(get(hObject,'String')) returns contents of dt_sim as a double


% --- Executes during object creation, after setting all properties.
function dt_sim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dt_sim (see GCBO)
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



function Nx_Callback(hObject, eventdata, handles)
% hObject    handle to Nx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nx as text
%        str2double(get(hObject,'String')) returns contents of Nx as a double


% --- Executes during object creation, after setting all properties.
function Nx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ny_Callback(hObject, eventdata, handles)
% hObject    handle to Ny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ny as text
%        str2double(get(hObject,'String')) returns contents of Ny as a double


% --- Executes during object creation, after setting all properties.
function Ny_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Nz_Callback(hObject, eventdata, handles)
% hObject    handle to Nz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nz as text
%        str2double(get(hObject,'String')) returns contents of Nz as a double


% --- Executes during object creation, after setting all properties.
function Nz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Nt_Callback(hObject, eventdata, handles)
% hObject    handle to Nt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nt as text
%        str2double(get(hObject,'String')) returns contents of Nt as a double


% --- Executes during object creation, after setting all properties.
function Nt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in make_grid.
function make_grid_Callback(hObject, eventdata, handles)
Nx = str2num(get(handles.Nx,'String'));
Ny = str2num(get(handles.Ny,'String'));
Nz = str2num(get(handles.Nz,'String'));
Nt = str2num(get(handles.Nt,'String'));
dx = str2num(get(handles.dx,'String'))*1e-3;
dy = str2num(get(handles.dy,'String'))*1e-3;
dz = str2num(get(handles.dz,'String'))*1e-3;
dt = str2num(get(handles.dt_sim,'String'))*1e-6;
if Nz == 0
    if Ny == 0
        kWaveGrid(Nx,dx);
    else
        kgrid = kWaveGrid(Nx,dx,Ny,dy);
    end
else
    kgrid = kWaveGrid(Nx,dx,Ny,dy,Nz,dz);
end
kgrid.Nt = Nt;
    kgrid.dt = dt;
kgrid.t_array = 0:dt:Nt*dt;
assignin('base','kgrid',kgrid);
% hObject    handle to make_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in make_medium.
%Input 0 for full range of xN, yN

function make_medium_Callback(hObject, eventdata, handles)
kgrid = evalin('base','kgrid');
mode = handles.medium_menu.String{handles.medium_menu.Value};
switch mode
    case 'Homogenous'
        sound_speed = str2double(handles.medium_ss.String)*1e3;
        alpha_coeff = str2double(handles.medium_alpha_coeff.String);
        alpha_power = str2double(handles.medium_alpha_power.String);
        density = str2double(handles.medium_density.String)*1e3;
        medium.sound_speed = sound_speed; %m/s
        medium.alpha_coeff = alpha_coeff; %dB/(MHz^y cm)
        medium.alpha_power = alpha_power;
        medium.density = density;
    case 'Heterogenous'
        BG = ones(kgrid.Nx,kgrid.Ny);
        sound_speed = str2double(handles.medium_ss.String)*1e3*BG;
        alpha_coeff = str2double(handles.medium_alpha_coeff.String)*BG;
        alpha_power = str2double(handles.medium_alpha_power.String);
        density = str2double(handles.medium_density.String)*1e3*BG;
        medium.sound_speed = sound_speed; %m/s
        medium.alpha_coeff = alpha_coeff; %dB/(MHz^y cm)
        medium.alpha_power = alpha_power;
        medium.density = density;
        %Layer 1
        if handles.ly1_on.Value
            xr1 = str2double(handles.ly1_xN.String);
            yr1 = str2double(handles.ly1_yN.String);
            if xr1 == 0
                xr1 = kgrid.Nx;
            end
            if yr1 == 0
                yr1 = kgrid.Ny;
            end
            rng1 = ones(xr1,yr1);
            ss1 = str2double(handles.ly1_ss.String)*rng1.*1e3;
            den1 = str2double(handles.ly1_den.String)*rng1.*1e3;
            ac1 = str2double(handles.ly1_ac.String)*rng1;
            x01 = str2double(handles.ly1_x0.String);
            y01 = str2double(handles.ly1_y0.String);
            medium.sound_speed(x01:x01+xr1-1,y01:y01+yr1-1) = ss1;
            medium.density(x01:x01+xr1-1,y01:y01+yr1-1) = den1;
            medium.alpha_coeff(x01:x01+xr1-1,y01:y01+yr1-1) = ac1;
        end
        %Layer 2
        if handles.ly2_on.Value
            xr1 = str2double(handles.ly2_xN.String);
            yr1 = str2double(handles.ly2_yN.String);
            if xr1 == 0
                xr1 = kgrid.Nx;
            end
            if yr1 == 0
                yr1 = kgrid.Ny;
            end
            rng1 = ones(xr1,yr1);
            ss1 = str2double(handles.ly2_ss.String)*rng1.*1e3;
            den1 = str2double(handles.ly2_den.String)*rng1.*1e3;
            ac1 = str2double(handles.ly2_ac.String)*rng1;
            x01 = str2double(handles.ly2_x0.String);
            y01 = str2double(handles.ly2_y0.String);
            medium.sound_speed(x01:x01+xr1-1,y01:y01+yr1-1) = ss1;
            medium.density(x01:x01+xr1-1,y01:y01+yr1-1) = den1;
            medium.alpha_coeff(x01:x01+xr1-1,y01:y01+yr1-1) = ac1;
        end
        %Layer 3
        if handles.ly3_on.Value
            xr1 = str2double(handles.ly3_xN.String);
            yr1 = str2double(handles.ly3_yN.String);
            if xr1 == 0
                xr1 = kgrid.Nx;
            end
            if yr1 == 0
                yr1 = kgrid.Ny;
            end
            rng1 = ones(xr1,yr1);
            ss1 = str2double(handles.ly3_ss.String)*rng1.*1e3;
            den1 = str2double(handles.ly3_den.String)*rng1.*1e3;
            ac1 = str2double(handles.ly3_ac.String)*rng1;
            x01 = str2double(handles.ly3_x0.String);
            y01 = str2double(handles.ly3_y0.String);
            medium.sound_speed(x01:x01+xr1-1,y01:y01+yr1-1) = ss1;
            medium.density(x01:x01+xr1-1,y01:y01+yr1-1) = den1;
            medium.alpha_coeff(x01:x01+xr1-1,y01:y01+yr1-1) = ac1;
        end
end

% f_max = medium.sound_speed/2/kgrid.dx; %maximum frequency allowed

assignin('base','medium',medium);
% hObject    handle to make_medium (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in make_source.
function make_source_Callback(hObject, eventdata, handles)
kgrid = evalin('base','kgrid');
if handles.source_3d.Value

    %get waveform and time grid
    type = handles.source_type1.String{handles.source_type1.Value};
    dt = kgrid.dt;
    t = kgrid.t_array*1e6;
    sampling_freq = 1/kgrid.dt;
    delay = str2num(get(handles.trans_delay,'String'));
    mag = str2num(get(handles.trans_mag,'String'));
    t0 = find(t >= delay,1);
    %get transducer characterstics
    steering_angle = str2double(get(handles.trans_steer,'String'));
    element_spacing_x = str2double(get(handles.trans_x_kerf,'String'))/1e3;
    element_spacing_y = str2double(get(handles.trans_y_kerf,'String'))/1e3;
    burst_freq = str2double(get(handles.trans_freq,'String'))*1e6;
    burst_cycles = str2double(get(handles.trans_cycles,'String'));
    num_elements = str2double(get(handles.trans_elements,'String'));

    %Create element Mask
    az_elements = str2double(handles.az_ele.String);
    ela_elements = str2double(handles.ela_ele.String);
    if num_elements ~= az_elements*ela_elements
        errordlg('Number of transducer elments must equal azi elements * elavational elements')
    end

    z_vec = kgrid.z_vec-min(kgrid.z_vec);
    y_vec = kgrid.y_vec-min(kgrid.y_vec);
    %%%%%%%%%%%%%%%% Create Basic Element Index

    y_kerf = str2double(handles.trans_x_kerf.String);
    z_kerf = str2double(handles.trans_y_kerf.String);

    y_loc = round(linspace(-az_elements*(y_kerf+1)/2,az_elements*(y_kerf+1)/2,az_elements));

    z_loc = round(linspace(-ela_elements*(z_kerf+1)/2,ela_elements*(z_kerf+1)/2,ela_elements));

    x0 = str2double(handles.trans_x0.String);
    y0 = str2double(handles.trans_y0.String);
    z0 = str2double(handles.trans_z0.String);

    y_element_index = y_loc+y0;
    x_element_index = x0;
    z_element_index = z0+z_loc;



    %%%%%%%%%%%%%%%%%%%%% Modify for element widths

    y_width = str2double(handles.trans_width.String);
    z_width = str2double(handles.trans_length.String);
    y_elements_width = az_elements*y_width;

    if y_width > 1
        y_space = str2double(handles.az_space.String);
        if mod(y_width,2) == 0
            y_element_addition = unique(sort([0, y_space:y_space:y_width*y_space/2, (y_space:y_space:y_width*y_space/2-y_space)*-1]));
        else

            y_element_addition = unique(sort([0, y_space:y_space:y_width*y_space/2, (y_space:y_space:y_width*y_space/2)*-1]));

        end
        for i = 1:az_elements
            for j = 1:y_width
                new_y_element_index((i-1)*y_width+j) = y_element_index(i)+y_element_addition(j);
            end
        end
    else
        new_y_element_index = y_element_index;
    end

    if z_width > 1
        z_space = str2double(handles.el_space.String);
        if mod(z_width,2) == 0
            z_element_addition = unique(sort([0, z_space:z_space:z_width*z_space/2, (z_space:z_space:z_width*z_space/2-z_space)*-1]));
        else
            z_element_addition = unique(sort([0, z_space:z_space:z_width*z_space/2, (z_space:z_space:z_width*z_space/2)*-1]));
        end
        for i = 1:ela_elements
            for j = 1:z_width
                new_z_element_index((i-1)*z_width+j) = z_element_index(i)+z_element_addition(j);
            end
        end
    else
        new_z_element_index = z_element_index;
    end

    %%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%  Finish Building Source Mask
    if ~handles.array_2d.Value
        all_element_index = zeros([kgrid.Nx,kgrid.Ny,kgrid.Nz]);
        for i = 1:kgrid.Nx
            for j = 1:kgrid.Ny
                for k = 1:kgrid.Nz
                    if ismember(i,x_element_index) && ismember(j,new_y_element_index) && ismember(k,new_z_element_index)
                        all_element_index(i,j,k) = 1;
                        %   else
                        %    all_element_index(i,j,k) = 0;
                    end
                end

            end

            multiWaitbar('Building Mask',i/kgrid.Nx);
        end
    else
        x_element_index = str2num(handles.curve_el.String);
        x_element_index = x_element_index+x0;
        if length(x_element_index) ~= ela_elements
            errordlg('Number of points in Elevational Dispalcement must equal number of elevational elements');
        end
        all_element_index = zeros([kgrid.Nx,kgrid.Ny,kgrid.Nz]);
        for i = 1:kgrid.Nx
            for j = 1:kgrid.Ny
                for k = 1:kgrid.Nz
                    if ismember(i,x_element_index) && ismember(j,new_y_element_index) && ismember(k,new_z_element_index)
                        z_spot = find(new_z_element_index == k);
                        z_spot = ceil(z_spot/z_width);
                        x_temp = x_element_index(z_spot);
                        if ismember(i,x_temp)
                            all_element_index(i,j,k) = 1;
                        end
                        % else
                        %        all_element_index(i,j,k) = 0;
                    end
                end

            end

            multiWaitbar('Building Mask',i/kgrid.Nx);
        end
    end
    %     multiWaitbar('CLOSEALL');
    %    voxelPlot(double(all_element_index))

    %%%%%%%%%%%%%%%%%% CREATE BURST OFFSETS
    % use geometric beam forming to calculate the tone burst offsets for each
    % transducer element based on the element index
    medium = evalin('base','medium');
    ss = str2double(handles.medium_ss.String)*1000;
    % element_space = kgrid.dy * (-(az_elements - 1)/2:(az_elements - 1)/2);
    element_space = kgrid.y_vec(y_element_index);
    % element_space = repmat(element_space,16,2);
    %             burst_offset = t0 + element_space * sin(steering_angle * pi/180)...
    %                 / (medium.sound_speed * kgrid.dt);


    wave_type = handles.wave_type.String{handles.wave_type.Value};
    switch wave_type
        case {'Plane','STA'}

%             burst_offset = t0 + element_space * sin(steering_angle * pi/180)/...
%                 (ss * kgrid.dt);
%             if handles.array_2d.Value
%                 a_center = max(x_element_index)-min(x_element_index);
%             else
%                 burst_offset  = repmat(burst_offset,ela_elements,1);
%             end
            fpoint = str2num(handles.trans_focus.String); %voxels
            fpoint(3) = 10000;
        case 'Focus'
            fpoint = str2num(handles.trans_focus.String); %voxels
    end
            x_point = fpoint(1);
            y_point = fpoint(2);
            z_point = fpoint(3);
            if handles.array_2d.Value

                rep_y = num_elements/az_elements;
                rep_z = num_elements/ela_elements;
                for i = 1:length(x_element_index)
                    x_element_index2((i-1)*rep_z+1:i*rep_z) = ones(1,rep_z).*x_element_index(i);
                    z_element_index2((i-1)*rep_z+1:i*rep_z) = ones(1,rep_z).*z_element_index(i);
                end
                %                                 x_element_index2 = repmat(x_element_index,rep_z);
                %                                 x_element_index2 = x_element_index2(1,:);
                y_element_index2 = repmat(y_element_index,rep_y);
                y_element_index2 = y_element_index2(1,:);
                %                                 z_element_index2 = repmat(z_element_index,rep_z);
                %                                 z_element_index2 = z_element_index2(1,:);

                %                                 for i = 1:length(x_element_index)

                all_element_loc = ones(3,num_elements);
                all_element_loc(2,:) = y_element_index2;
                all_element_loc(1,:) = x_element_index2;
                all_element_loc(3,:) = z_element_index2;
                if strcmp(wave_type,'Plane')
                    aloc2 = all_element_loc(2,:)-fpoint(2);
                    aloc2 = round(aloc2*sin(steering_angle * pi/180));
                end

                %                                 all_element_loc(1,:) = sort(all_element_loc(1,:));
                %                                 all_element_loc(3,:) = sort(all_element_loc(3,:));

                if x_point == x0
                    for i = 1:num_elements
                        foc_dist(i) = (x_point-all_element_loc(1,i).*kgrid.dx);
                    end
                else
                    for i = 1:num_elements
                        if strcmp(wave_type,'Plane')
                            foc_dist(i) = sqrt(((x_point-all_element_loc(1,i)).*kgrid.dx).^2+((y_point-all_element_loc(2,i)+aloc2(i)).*kgrid.dy).^2+((z_point-all_element_loc(3,i)).*kgrid.dz).^2);
                        else
                            foc_dist(i) = sqrt(((x_point-all_element_loc(1,i)).*kgrid.dx).^2+((y_point-all_element_loc(2,i)).*kgrid.dy).^2+((z_point-all_element_loc(3,i)).*kgrid.dz).^2);
                        end
                    end
                end
            else
                if num_elements == az_elements
                    all_element_loc = ones(3,num_elements);
                    all_element_loc(2,:) = y_element_index;
                    all_element_loc(1,:) = ones(1,num_elements).*x_element_index;
                    all_element_loc(3,:) = ones(1,num_elements).*z_element_index;
                else
                    rep_y = num_elements/az_elements;
                    rep_z = num_elements/ela_elements;
                    y_element_index2 = repmat(y_element_index,rep_y);
                    y_element_index2 = y_element_index2(1,:);
                    z_element_index2 = repmat(z_element_index,rep_z);
                    z_element_index2 = z_element_index2(1,:);
                    all_element_loc = ones(3,num_elements);
                    all_element_loc(2,:) = y_element_index2;
                    all_element_loc(1,:) = ones(1,num_elements).*x_element_index;
                    all_element_loc(3,:) = sort(z_element_index2);
                end

                for i = 1:num_elements
                    foc_dist(i) = sqrt(((x_point-all_element_loc(1,i)).*kgrid.dx).^2+((y_point-all_element_loc(2,i)).*kgrid.dy).^2+((z_point-all_element_loc(3,i)).*kgrid.dz).^2);
                end
            end
    



    %                         e_ind = element_index*(str2num(get(handles.trans_x_kerf,'String'))+elwidth);
    %                         w = fpoint(2)-y0; %depth in mm/10
    %                         h = (e_ind-(fpoint(1)-x0))*-1; %Height in mm/10
    %                         hyp = sqrt(h.^2+w^2); %hypoteneuse in mm/10
    %                         S = hyp./(ss*kgrid.dt*1e4); %hypoteneuse in samples
    if handles.usemedia.Value
        % h = round(h);

        for i = 1:num_elements
            d = all_element_loc(1,i);
            theta(1,i) = atan((fpoint(1)-d)/(all_element_loc(2,i)-fpoint(2)));
            theta(2,i) = atan((fpoint(1)-d)/(all_element_loc(3,i)-fpoint(3)));
        end

        for i = 1:fpoint(1)-x0
            for j = 1:num_elements
                sslocy(i,j) = all_element_loc(2,j)-round(i/tan(theta(1,j)));
                sslocz(i,j) = all_element_loc(3,j)-round(i/tan(theta(2,j)));
            end
        end
        for i = 1:fpoint(1)-x0
            for j = 1:num_elements
                ss2(i,j) = medium.sound_speed(i+x0,sslocy(i,j),sslocz(i,j));
            end
        end
        ss3 = mean(ss2,1);
        %   hyp2 = round(h);
        burst_offset = foc_dist./ss3./kgrid.dt; %S;
        burst_offset = abs(burst_offset-min(burst_offset));
    else
        if length(medium.sound_speed) > 1
            medium.sound_speed = medium.sound_speed(1);
        end
        burst_offset = foc_dist/medium.sound_speed(1)/kgrid.dt;
        burst_offset = abs(burst_offset-min(burst_offset));
        burst_offset = fliplr(burst_offset);
    end
    %                        burst_offset = S;
    %                         burst_offset = foc_dist/medium.sound_speed(1)/kgrid.dt; %added for 3D. Havent included heterogenous yet
    %                         burst_offset = abs(burst_offset-min(burst_offset));
    %                 for i = 1:length(h)
    %                     hyp = sqrt(h.^2+w^2);
    %                 end
    %                     errordlg('You are an ass hole')

    % Plot Burst Offset
    axes(handles.axes2); plot(burst_offset*kgrid.dt*1e6*-1)
    xlabel('Element Number')
    ylabel('Delay \mus')

    % Create tone burst signals
    slength = kgrid.Nt-101; %change the subtractor for different length signals
    source.p = mag * toneBurst(sampling_freq, burst_freq, burst_cycles,...
        'SignalOffset',burst_offset,'SignalLength',slength);

    %%%%%%%%%%% ADJUST FOR ELEMENT WIDTH %%%%%%%%%%%
    source.p2 = zeros(num_elements*y_width,size(source.p,2));
    for i = 1:num_elements
        for j = 1:y_width%*z_width
            %for k = 1:z_width
            source.p2((i-1)*y_width+j,:) = source.p(i,:);
            %  end
        end
        multiWaitbar('Creating source structure in y',i/num_elements);
    end
    if handles.array_2d.Value
        num_elements2 = num_elements.*y_width;
        num_elements3 = num_elements./z_width;
        for i = 1:num_elements2
            %                         for i  = 1:num_elements3
            %  for k = 1:az_elements
            for j = 1:z_width%*z_width
                %for k = 1:z_width
                %      source.p3((i-1)*az_elements+k,:) = source.p2(i,:);
                source.p3((i-1)*z_width+j,:) = source.p2(i,:);
                % source.p3((j-1)*num_elements+i,:) = source.p2(i,:);
                %  end
            end
            multiWaitbar('Creating source structure in z',i/num_elements2);
        end
        % new_z = 1:az_elements;
        if y_width > 1
            new_z = 1:az_elements*y_width;
            for i = 1:az_elements
                for j = 1:y_width
                    old_z(((i-1)*y_width+j)) = (i-1)*y_width*z_width+j;
                end
            end
        else
            new_z = 1:(az_elements); %adjust to az_elements*y_width for width transducers?
            old_z = new_z.*z_width*y_width-(z_width*y_width-1);
        end


        for i = 1:length(x_element_index)
            for j = 1:z_width
                %                               for k = 1:y_width
                % source.p4(1+(j-1)*az_elements*i+k,:) = source.p3()
                % source.p4{i}(1+(j-1)*az_elements+k,:) = source.p3(i*az_elements*z_width+,:)
                source.p4{i}(new_z+(j-1)*az_elements*y_width,:) = source.p3(old_z+(j-1)*y_width+(i-1)*(az_elements*z_width*y_width),:);
                %                                   source.p4{i}(new_z+(j-1)*az_elements,:) = source.p3(old_z+(j-1)+(i-1)*(az_elements*z_width),:);
                %                               end
            end
        end
        for i = 1:length(source.p4)
            source.p3(1+(i-1)*size(source.p4{1},1):i*size(source.p4{1},1),:) = source.p4{i};
        end

        source.p2 = source.p3;
    else
        num_y = size(source.p2,1);
        if z_width > 1
            for i = 1:z_width
                source.p3((i-1)*num_y+1:i*num_y,:) = source.p2(:,:);
                multiWaitbar('Creating source structure i z',i/z_width);
            end

            source.p2 = source.p3;
        end
    end
    %
    source.p = source.p2;


    if strcmp('Focus',handles.wave_type.String{handles.wave_type.Value})
        source.p = fliplr(source.p);
        for i = 1:size(source.p,1)
            s(i) = find(source.p(i,:),1)-1;
        end
        minS = min(s);
        source.p = circshift(source.p,-minS,2);
    end

    source_1 = zeros(kgrid.Nx,kgrid.Ny,kgrid.Nz);

    %Create mask
    %
    source.p_mask = all_element_index;
    if ~handles.had_box.Value
        active_elements = str2num(get(handles.trans_active_elements,'String'));
        if num_elements == az_elements
            element_active = y_element_index(str2num(get(handles.trans_active_elements,'String')));
        else
            element_active = y_element_index2(active_elements);
        end
        if length(element_active) < num_elements
            active_elements = active_elements.*y_kerf;
            %                           mask_size = size(all_element_index);
            %                           all_element_index = reshape(all_element_index,[numel(all_element_index),1]);
            %                           active_elements = find(all_element_index == 1);
            %  element_active2 = active_elements(element_active);
            if y_width > 1
                if mod(y_width,2) == 0
                    y_element_addition = unique(sort([0, y_space:y_space:y_width*y_space/2, (y_space:y_space:y_width*y_space/2-y_space)*-1]));
                else
                    y_element_addition = unique(sort([0, y_space:y_space:y_width*y_space/2, (y_space:y_space:y_width*y_space/2)*-1]));
                end
                for i = 1:length(element_active)
                    for j = 1:y_width
                        e_active4((i-1)*y_width*z_width+j) = element_active(i)+y_element_addition(j);
                        active_elements2((i-1)*y_width*z_width+j) = active_elements(i)+y_element_addition(j);
                        % element_active3((i-1)*y_width*z_width+j) = element_active2(i)+y_element_addition(j)*mask_size(1);
                        %   element_active3((i-1)*y_width*z_width+j) = element_active2(i)+y_element_addition(j);
                    end
                end
                %  element_active2 = element_active3;
                element_active = e_active4;
                active_elements = active_elements2;
            end
            %                           all_element_index = zeros(size(all_element_index));
            %                           all_element_index(element_active2) = 1;
            %                           source.p_mask = reshape(all_element_index,[mask_size]);
            all_element_index = zeros([kgrid.Nx,kgrid.Ny,kgrid.Nz]);
            for i = 1:kgrid.Nx
                for j = 1:kgrid.Ny
                    for k = 1:kgrid.Nz
                        if ismember(i,x_element_index) && ismember(j,element_active) && ismember(k,z_element_index)
                            all_element_index(i,j,k) = 1;
                        else
                            all_element_index(i,j,k) = 0;
                        end
                    end

                end

                multiWaitbar('Restructuring Mask for Active Elements',i/kgrid.Nx);
            end

            source.p_mask = all_element_index;

            source.p2 = source.p(active_elements-1,:);
            clear source.p;
            source.p = source.p2;
            clear source.p2;
        end

    else %Hadamard
        %  if handles.wave_type.Value == 3
        sta = str2num(handles.sta.String);
        groups = str2double(handles.new_foci.String);
        sta2 = ones(size(source.p,1),1);
%         for i = 1:ela_elements
%             for j = 1:az_elements
%                 source.p2((1+(j-1)*az_elements*y_width*z_width+(i-1)*ela_elements):(j*az_elements*y_width*z_width+(i-1)*ela_elements)) = ...
%                     source.p((1+(j-1)*az_elements*y_width*z_width+(i-1)*ela_elements):(j*az_elements*y_width*z_width+(i-1)*ela_elements)).*sta(j);
%             end
%         end
        for k = 1:az_elements
%             for i = 1:y_width
%                 for j = 1:z_width
                    sta3(1+(k-1)*y_width*z_width:k*y_width*z_width) = sta(k);
%                 end
%             end
        end
        for i = 1:ela_elements
            sta4(1+(length(sta3)*(i-1)):length(sta3)*i) = sta3;
        end
        for i = 1:size(source.p,1)
            source.p2(i,:) = source.p(i,:).*sta4(i);
        end
                    %         for i = 1:ela_elements
                    %             for j = 1:z_width
                    %                 for k = 1:az_elements
                    %                     source.p2((k+(i-1)*ela_elements):(k+(i-1)*ela_elements+y_width)) = source.p(k+(i-1)*ela_elements)*sta(k);
%                 end
%             end
% %                  source.p2((1+(k-1)*az_elements+(i-1)*ela_elements+(j-1)*z_width):
% %             source.p2((i-1)*groups+1,:) = source.p((i-1)*groups+1,:).*sta(i);
%         end
        source.p = source.p2;
    end
    if handles.negfocus.Value
        for i = 1:size(source.p,1)
            t3(i) = find(source.p(i,:),1);
        end
        t2 = abs(max(t3)-t3);
        %Create new plane wave source.p
        steering_angle = 0;
        burst_offset = t0 + element_space * sin(steering_angle * pi/180)/...
            (ss * kgrid.dt);
        source.p = mag * toneBurst(sampling_freq, burst_freq, burst_cycles,...
            'SignalOffset',burst_offset(1:size(source.p2,1)),'SignalLength',slength);
        %shift source.p to the inverse of t0
        for i = 1:size(source.p,1)
            %                      if ismember(i,element_active2)
            source.p(i,:) = circshift(source.p(i,:),t2(i));
            %                      else
            %                          source.p(i,:) = zeros(1,size(source.p,2));
            %                      end
        end
    end
    multiWaitbar('CLOSEALL');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else %%%%%%%%%%%%%%%%%    2D     %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get waveform and time grid
    type = handles.source_type1.String{handles.source_type1.Value};
    dt = kgrid.dt;
    t = kgrid.t_array*1e6;
    sampling_freq = 1/kgrid.dt;
    delay = str2num(get(handles.trans_delay,'String'));
    mag = str2num(get(handles.trans_mag,'String'));
    t0 = find(t >= delay,1);
    %get transducer characterstics

    steering_angle = str2double(get(handles.trans_steer,'String'));

    element_spacing = str2double(get(handles.trans_x_kerf,'String'))/1e3;
    burst_freq = str2double(get(handles.trans_freq,'String'))*1e6;
    burst_cycles = str2double(get(handles.trans_cycles,'String'));
    num_elements = str2double(get(handles.trans_elements,'String'));
    x_vec = kgrid.x_vec-min(kgrid.x_vec);
    %             for i = 1:num_elements
    %                 [~,element_index(i)] = min(abs(x_vec - i*element_spacing));
    %             end
    element_index = -(num_elements -1)/2:(num_elements -1)/2;

    % use geometric beam forming to calculate the tone burst offsets for each
    % transducer element based on the element index
    medium = evalin('base','medium');
    ss = str2double(handles.medium_ss.String)*1000;
    x_kerf = str2double(handles.trans_x_kerf.String);
    element_space = kgrid.dx * (-(num_elements - 1)/2:(num_elements - 1)/2) * x_kerf;
    %             burst_offset = t0 + element_space * sin(steering_angle * pi/180)...
    %                 / (medium.sound_speed * kgrid.dt);

    %THIS IS WHERE FOCUS AND PLANE DELAYS ARE CREATED
    width = str2num(get(handles.trans_width,'String')); %voxels
    elwidth = str2num(get(handles.trans_length,'String')); %voxels
    wave_type = handles.wave_type.String{handles.wave_type.Value};
    switch wave_type
        case 'Plane'


            for j = 1:length(element_space)
                %  if t(i) <= 0
                delays(j) = tand(steering_angle).*element_space(j)/(ss*kgrid.dt);
                % else
                %    delays(i,j) = tand(t(i)).*abs(xloc3(end)-xloc3(j));
                % end
            end


            delays = delays-min(delays);

            %                             for i = 1:length(delays)
            %                                 delays(i) = 1 + abs(delays(i)-max(delays))/(ss*kgrid.dt);
            %                             end

            burst_offset = t0 + element_space * sin(steering_angle*pi/180)/...
                (ss * kgrid.dt);
            burst_offset = burst_offset - min(burst_offset);
        case 'Focus' %Focus
            x0 = str2double(handles.trans_x0.String); %voxels
            y0 = str2double(handles.trans_y0.String); %voxels
            fpoint = str2num(handles.trans_focus.String); %voxels
            %                 for i = 1:length(element_index)
            %                     steering_angle = atand((fpoint(2)-y0)/(element_index(i) + x0-fpoint(1)));
            %                     burst1 = t0 + element_space * sin(steering_angle * pi/180)/...
            %                         (ss * kgrid.dt);
            %                     burst_offset(i) = burst1(i);
            %                 end
            e_ind = element_index*(str2num(get(handles.trans_x_kerf,'String'))+width);
            w = fpoint(2)-y0; %depth in mm/10
            h = (e_ind-(fpoint(1)-x0))*-1; %Height in mm/10
            hyp = sqrt((h.*kgrid.dx).^2+(w.*kgrid.dy)^2); %hypoteneuse in mm/10
            S = hyp./(ss*kgrid.dt/1e4); %hypoteneuse in samples
            if handles.usemedia.Value
                h = round(h);

                for i = 1:length(h)
                    theta(i) = atan(w/h(i));
                end

                for i = 1:w
                    for j = 1:length(h)
                        sslocx(i,j) = h(j)-round(i/tan(theta(j)));
                    end
                end
                for i = 1:w
                    for j = 1:length(h)
                        ss2(i,j) = medium.sound_speed(sslocx(i,j)+x0,i+y0);
                    end
                end
                ss3 = mean(ss2,1);
                %   hyp2 = round(h);
                burst_offset = hyp./ss3./kgrid.dt; %S;
                burst_offset = abs(burst_offset-min(burst_offset))+1;
                burst_offset = fliplr(burst_offset);
                %                     hyp2 = hyp;
                %                                 for i = 1:length(hyp2)
                %                                     xloc{i} = round(linspace(x0-h(i),fpoint(1),hyp2(i)));
                %                                     yloc{i} = round(linspace(y0,fpoint(2),hyp2(i)));
                %                                 end
                %                                 for i = 1:length(hyp2)
                %                                     for j = 1:hyp2(i)
                %                                         pat{i}(j) = medium.sound_speed(round(xloc{i}(j)*(0.0001/kgrid.dx)),round(yloc{i}(j)*(0.0001/kgrid.dy)));
                %                                     end
                %                                 end
                %                                 for i = 1:length(pat)
                %                                     ss2(i) = sum(pat{i})/hyp(i);
                %                                 end
                %                                 S = hyp./(ss2*kgrid.dt)/1e4;
            else
                if length(medium.sound_speed) > 1
                    medium.sound_speed = medium.sound_speed(1);
                end
                burst_offset = hyp/medium.sound_speed/kgrid.dt; %S;
                burst_offset = abs(burst_offset-min(burst_offset))+1;
                %                                 burst_offset = fliplr(burst_offset);
            end

            %                 for i = 1:length(h)
            %                     hyp = sqrt(h.^2+w^2);
            %                 end
            %                     errordlg('You are an ass hole')
            %                         case 'STA'
            %                              burst_offset = t0 + element_space * sin(steering_angle * pi/180)/...
            %                             (ss * kgrid.dt);
    end
    % Plot Burst Offset
    axes(handles.axes2); plot(burst_offset*kgrid.dt*1e6*-1)
    xlabel('Element Number')
    ylabel('Delay \mus')

    % Create tone burst signals
    slength = kgrid.Nt-101; %change the subtractor for different length signals
    source.p = mag * toneBurst(sampling_freq, burst_freq, burst_cycles,...
        'SignalOffset',burst_offset,'SignalLength',slength);
    switch wave_type
        case 'Focus'
            source.p = fliplr(source.p);
            %                 if handles.gauss_apod.Value
            %                     H = gausswin(size(source.p,1));
            %                     shiftpt = find(source.p,1,'last');
            for i = 1:size(source.p,1)
                s(i) = find(source.p(i,:),1)-1;
            end
            minS = min(s);
            source.p = circshift(source.p,-minS,2);

            %                             %


    end

    source_1 = zeros(kgrid.Nx,kgrid.Ny);

    %Create mask
    x0 = str2double(get(handles.trans_x0,'String'));
    element_index2 = round(x0 + element_index * (element_spacing + 1/1000)* 1e3);
    y0 = str2double(get(handles.trans_y0,'String'));
    source.p_mask = zeros(kgrid.Nx,kgrid.Ny);
    element_active = str2num(get(handles.trans_active_elements,'String'));
    %             source.p_mask(element_index2,y0) = 1;

    %%%%%%%% Adjust for element width and active elements
    %%%%%%%% %%%%%%%%%%%
    active = element_index2(element_active);
    %%%%
    y_space = str2double(handles.az_space.String);
    if mod(width,2) == 0
        w_add = -width*y_space/2+y_space:y_space:width*y_space/2;
    else
        if mod(y_space,2) == 0
            w_add = floor(width*y_space/2)*-1+1:y_space:floor(width*y_space/2)-1;
        else
            w_add = floor(width*y_space/2)*-1:y_space:floor(width*y_space/2);
        end
    end
    for i = 1:length(active)
        for j = 1:width
            active2((i-1)*width+j) = active(i)+w_add(j);
        end
    end


    %                     w3 = -w2:w2;
    %                     w4 = 1:width;

    for i = 1:kgrid.Nx
        for j = 1:kgrid.Ny
            if ismember(i,active2) && ismember(j,y0)
                source.p_mask(i,j) = 1;
            end
        end
    end
    %                   
    source.p = source.p(element_active,:);
    for i = 1:length(active)
        for j = 1:width
            source.p2((i-1)*width+j,:) = source.p(element_active(i),:);
        end
    end
    clear source.p;
    source.p = source.p2;
    clear source.p2;
    if handles.had_box.Value
        % if handles.wave_type.Value == 3
        sta = str2num(handles.sta.String);
        for i = 1:length(sta)
            for j = 1:width
                source.p((i-1)*width+j,:) = source.p((i-1)*width+j,:).*sta(i);
            end
        end
    end
    if handles.negfocus.Value
        for i = 1:size(source.p,1)
            t3(i) = find(source.p(i,:),1);
        end
        t2 = abs(max(t3)-t3);
        %Create new plane wave source.p
        steering_angle = 0;
        burst_offset = t0 + element_space * sin(steering_angle * pi/180)/...
            (ss * kgrid.dt);
        source.p = mag * toneBurst(sampling_freq, burst_freq, burst_cycles,...
            'SignalOffset',burst_offset(1:size(source.p2,1)),'SignalLength',slength);
        %shift source.p to the inverse of t0
        for i = 1:size(source.p,1)
            %                      if ismember(i,element_active2)
            source.p(i,:) = circshift(source.p(i,:),t2(i));
            %                      else
            %                          source.p(i,:) = zeros(1,size(source.p,2));
            %                      end
        end
    end


end
n1 = size(source.p,1);
n2 = length(find(source.p_mask));
set(handles.source_n,'String',num2str(n1));
set(handles.mask_n,'String',num2str(n2));
if n2 == n1
    set(handles.mask_n,'ForegroundColor','k')
else
    set(handles.mask_n,'ForegroundColor','r')
end


clear source.p0
assignin('base','source',source);
% hObject    handle to make_source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in make_sensor.
function make_sensor_Callback(hObject, eventdata, handles)
% rad = str2num(get(handles.sensor_rad,'String'));
% pts = str2num(get(handles.sensor_pts,'String'));

if ismember('kgrid',evalin('base','who'))
%     mode_string = get(handles.sensor_mode,'String');
%     mode_val = get(handles.sensor_mode,'Value');
%     mode = mode_string{mode_val};
       rec_p = {'p'};
    if handles.sensor_velocity.Value
        rec_u = {'u'};
    else
        rec_u = {};
    end
    if handles.sensor_final.Value
        if handles.sensor_velocity.Value
        rec_f = {'p_final','u_final'};
        else 
            rec_f = {'p_final'};
        end
    else 
        rec_f = {}; 
    end
    if handles.sensor_p_max.Value
        rec_p_max = {'p_max'};
    else
        rec_p_max = {};
    end
    if handles.sensor_p_max.Value
        rec_p_max_all = {'p_max_all'};
    else
        rec_p_max_all = {};
    end
    if handles.sensor_intensity.Value
        rec_intensity = {'I'};
    else
        rec_intensity = {};
    end
    sensor.record = [rec_p, rec_u, rec_f, rec_p_max, rec_p_max_all, rec_intensity];
%%%%%%%%%%%%%%%%%%%%%%%%%% Imported from ksim_sensor %%%%%%%%%%%%%%%%%%%%%%
z = str2double(handles.s_plane.String);
p = handles.p_menu.String{handles.p_menu.Value};
kgrid = evalin('base','kgrid');
if ~handles.all_sensor.Value
    if handles.box3d.Value
        sizes = size(kgrid.x);
        switch p
            case 'XY'
                plane = sizes(1).*sizes(2);
                x1 = kgrid.x(:,:,z);
                y1 = kgrid.y(:,:,z);
                z1 = kgrid.z(:,:,z);
                x = reshape(x1,[1,plane]);
                y = reshape(y1,[1,plane]);
                z = reshape(z1,[1,plane]);
            case 'XZ'
                plane = sizes(1).*sizes(3);
                x1 = kgrid.x(:,z,:);
                y1 = kgrid.y(:,z,:);
                z1 = kgrid.z(:,z,:);
                x = reshape(x1,[1,plane]);
                y = reshape(y1,[1,plane]);
                z = reshape(z1,[1,plane]);
            case 'YZ'
                plane = sizes(2).*sizes(3);
                x1 = kgrid.x(z,:,:);
                y1 = kgrid.y(z,:,:);
                z1 = kgrid.z(z,:,:);
                x = reshape(x1,[1,plane]);
                y = reshape(y1,[1,plane]);
                z = reshape(z1,[1,plane]);
        end

    else
        x = linspace(1,kgrid.Nx,kgrid.Nx);
        for i = 1:length(y)
            y1(:,i) = y(i)*ones(1,kgrid.Ny);
            x1(:,i) = x;
        end
        y2 = reshape(y1,[1,numel(y1)]);
        x2 = reshape(x1,[1,numel(x1)]);
        y = y2;
        x = x2;
    end

else %All
    x = reshape(kgrid.x,[1,numel(kgrid.x)]);
    y = reshape(kgrid.y,[1,numel(kgrid.y)]);
    if handles.box3d.Value
        z = reshape(kgrid.z,[1,numel(kgrid.z)]);
    end
end
if kgrid.Nz == 0
    clear sensor.mask
    sensor.mask(1,:) = x;
    sensor.mask(2,:) = y;
else
    clear sensor.mask
    sensor.mask(1,:) = x;
    sensor.mask(2,:) = y;
    sensor.mask(3,:) = z;

    %sensor.mask(:,y,z) = x;
end
if handles.box3d.Value
    sensor.dims = 3;
else
    sensor.dims = 2;
end
if handles.all_sensor.Value
    sensor.all = 1;
else
    sensor.all = 0;
    plane = handles.p_menu.String{handles.p_menu.Value};
    sensor.plane = plane;
end
    assignin('base','sensor',sensor);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     
else
    errordlg('Need to create Kgrid before sensor')  
end
% hObject    handle to make_sensor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function sensor_rad_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensor_rad as text
%        str2double(get(hObject,'String')) returns contents of sensor_rad as a double


% --- Executes during object creation, after setting all properties.
function sensor_rad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensor_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sensor_pts_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_pts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensor_pts as text
%        str2double(get(hObject,'String')) returns contents of sensor_pts as a double


% --- Executes during object creation, after setting all properties.
function sensor_pts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensor_pts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in sensor_mode.
function sensor_mode_Callback(hObject, eventdata, handles)

% hObject    handle to sensor_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sensor_mode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sensor_mode


% --- Executes during object creation, after setting all properties.
function sensor_mode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensor_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sensor_xpos_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_xpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensor_xpos as text
%        str2double(get(hObject,'String')) returns contents of sensor_xpos as a double


% --- Executes during object creation, after setting all properties.
function sensor_xpos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensor_xpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sensor_ypos_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_ypos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensor_ypos as text
%        str2double(get(hObject,'String')) returns contents of sensor_ypos as a double


% --- Executes during object creation, after setting all properties.
function sensor_ypos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensor_ypos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in center_sensor.
function center_sensor_Callback(hObject, eventdata, handles)
if ismember('kgrid',evalin('base','who'))
    switch handles.sensor_mode.String{handles.sensor_mode.Value};
        case 'Cart_Circle'
              set(handles.sensor_xpos,'String',0);
    set(handles.sensor_ypos,'String',0);    
        case 'Circle'
    kgrid = evalin('base','kgrid');
    Ny = kgrid.Ny;
    Nx = kgrid.Nx;
    x = ceil(Nx/2);
    y = ceil(Ny/2);
    set(handles.sensor_xpos,'String',x);
    set(handles.sensor_ypos,'String',y);
    end
else
    errordlg('Need to create Kgrid before sensor')  
end
% hObject    handle to center_sensor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in reorder.
function reorder_Callback(hObject, eventdata, handles)
kgrid = evalin('base','kgrid');
sensor = evalin('base','sensor');
sensor_data = evalin('base','sensor_data');
sensor_data_reordered = reorderSensorData(kgrid, sensor, sensor_data);
assignin('base','sensor_data',sensor_data_reordered);
% hObject    handle to reorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function sensor_ypts_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_ypts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensor_ypts as text
%        str2double(get(hObject,'String')) returns contents of sensor_ypts as a double


% --- Executes during object creation, after setting all properties.
function sensor_ypts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensor_ypts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sensor_xpts_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_xpts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensor_xpts as text
%        str2double(get(hObject,'String')) returns contents of sensor_xpts as a double


% --- Executes during object creation, after setting all properties.
function sensor_xpts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensor_xpts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in source_menu1.
function source_menu1_Callback(hObject, eventdata, handles)
kgrid = evalin('base','kgrid');
switch handles.source_menu1.String{handles.source_menu1.Value}
    case 'none'
        set(handles.source_text1,'Visible','off')
        set(handles.source_text2,'Visible','off')
        set(handles.source_text3,'Visible','off')
        set(handles.source_text4,'Visible','off')
        set(handles.source_text13,'Visible','off')
        set(handles.source_text14,'Visible','off')
        set(handles.source_text15,'Visible','off')
        set(handles.source_text16,'Visible','off')
        set(handles.source_mag1,'Visible','off')
        set(handles.source_xpos1,'Visible','off')
        set(handles.source_ypos1,'Visible','off')
        set(handles.source_rad1,'Visible','off')
        set(handles.source_delay1,'Visible','off')
        set(handles.source_duration1,'Visible','off')
        set(handles.source_type1,'Visible','off')
        set(handles.source_period1,'Visible','off')
        set(handles.us_transducer,'Visible','off');
        set(handles.sourcegroup1,'Visible','off');
             set(handles.source_menu2,'Visible','on');
        set(handles.source_menu3,'Visible','on');
        set(handles.text22,'Visible','on');
        set(handles.text21,'Visible','on');
    case 'Disk'
        set(handles.source_text1,'Visible','on')
        set(handles.source_text2,'Visible','on')
        set(handles.source_text3,'Visible','on')
        set(handles.source_text4,'Visible','on')
        set(handles.source_text13,'Visible','off')
        set(handles.source_text14,'Visible','off')
        set(handles.source_text15,'Visible','off')
         set(handles.source_text16,'Visible','off')
        set(handles.source_mag1,'Visible','on')
        set(handles.source_xpos1,'Visible','on')
        set(handles.source_ypos1,'Visible','on')
        set(handles.source_rad1,'Visible','on')
        set(handles.source_delay1,'Visible','off')
        set(handles.source_duration1,'Visible','off')
        set(handles.source_type1,'Visible','off')
        set(handles.source_period1,'Visible','off')
        set(handles.us_transducer,'Visible','off');
        set(handles.sourcegroup1,'Visible','on');
             set(handles.source_menu2,'Visible','on');
        set(handles.source_menu3,'Visible','on');
        set(handles.text22,'Visible','on');
        set(handles.text21,'Visible','on');
        if kgrid.Nz > 0
            set(handles.source_zpos1,'Visible','on')
        else
            set(handles.source_zpos1,'Visible','off')
        end
    case 'Time-Varying Point'
        set(handles.source_text1,'Visible','on')
        set(handles.source_text2,'Visible','on')
        set(handles.source_text3,'Visible','on')
        set(handles.source_text4,'Visible','off')
        set(handles.source_text13,'Visible','on')
        set(handles.source_text14,'Visible','on')
        set(handles.source_text15,'Visible','on')
        set(handles.source_text16,'Visible','on')
        set(handles.source_mag1,'Visible','on')
        set(handles.source_xpos1,'Visible','on')
        set(handles.source_ypos1,'Visible','on')
        set(handles.source_rad1,'Visible','off')
        set(handles.source_delay1,'Visible','on')
        set(handles.source_duration1,'Visible','on')
        set(handles.source_type1,'Visible','on')
        set(handles.source_period1,'Visible','on')
        set(handles.us_transducer,'Visible','off');
        set(handles.sourcegroup1,'Visible','on');
        set(handles.source_menu2,'Visible','on');
        set(handles.source_menu3,'Visible','on');
        set(handles.text22,'Visible','on');
        set(handles.text21,'Visible','on');
        if kgrid.Nz > 0
            set(handles.source_zpos1,'Visible','on')
        else
            set(handles.source_zpos1,'Visible','off')
        end
        set(handles.trans_duration','Visible','on');
        set(handles.trans_period','Visible','on');
        set(handles.trans_mag','Visible','on');
        set(handles.trans_type,'Visible','on');
    case 'Transducer'
          set(handles.source_text1,'Visible','on')
        set(handles.source_text2,'Visible','on')
        set(handles.source_text3,'Visible','on')
        set(handles.source_text4,'Visible','on')
        set(handles.source_text13,'Visible','on')
        set(handles.source_text14,'Visible','on')
        set(handles.source_text15,'Visible','on')
        set(handles.source_text16,'Visible','on')
        set(handles.source_mag1,'Visible','on')
        set(handles.source_xpos1,'Visible','on')
        set(handles.source_ypos1,'Visible','on')
        set(handles.source_rad1,'Visible','on')
        set(handles.source_delay1,'Visible','on')
        set(handles.source_duration1,'Visible','on')
        set(handles.source_type1,'Visible','on')
        set(handles.source_period1,'Visible','on')
        set(handles.source_menu2,'Value',1);
        set(handles.source_menu3,'Value',1);
        set(handles.source_menu2,'Visible','off');
        set(handles.source_menu3,'Visible','off');
        set(handles.text22,'Visible','off');
        set(handles.text21,'Visible','off');
        set(handles.us_transducer,'Visible','on');
        set(handles.sourcegroup1,'Visible','off');
        set(handles.sourcegroup2,'Visible','off');
        set(handles.sourcegroup3,'Visible','off');
        set(handles.trans_arc_pos,'Visible','off');
        set(handles.trans_radius,'Visible','off');
        set(handles.trans_diameter','Visible','off');
        set(handles.trans_duration','Visible','off');
        set(handles.trans_period','Visible','off');
        set(handles.trans_mag','Visible','off');
         set(handles.trans_type,'Visible','off');
         set(handles.gauss_apod,'Visible','off');
         set(handles.negfocus,'Visible','off');
        
end
% hObject    handle to source_menu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns source_menu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from source_menu1


% --- Executes during object creation, after setting all properties.
function source_menu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_menu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in source_menu2.
function source_menu2_Callback(hObject, eventdata, handles)
kgrid = evalin('base','kgrid');
switch hObject.String{hObject.Value}
    case 'none'
        set(handles.source_text5,'Visible','off')
        set(handles.source_text6,'Visible','off')
        set(handles.source_text7,'Visible','off')
        set(handles.source_text8,'Visible','off')
        set(handles.source_mag2,'Visible','off')
        set(handles.source_xpos2,'Visible','off')
        set(handles.source_ypos2,'Visible','off')
        set(handles.source_rad2,'Visible','off')
        set(handles.sourcegroup2,'Visible','off');
    case 'Disk'
          set(handles.source_text5,'Visible','on')
        set(handles.source_text6,'Visible','on')
        set(handles.source_text7,'Visible','on')
        set(handles.source_text8,'Visible','on')
        set(handles.source_mag2,'Visible','on')
        set(handles.source_xpos2,'Visible','on')
        set(handles.source_ypos2,'Visible','on')
        set(handles.source_rad2,'Visible','on')
         set(handles.sourcegroup2,'Visible','on');
           if kgrid.Nz > 0
            set(handles.source_zpos2,'Visible','on')
        else
            set(handles.source_zpos2,'Visible','off')
        end
end
% hObject    handle to source_menu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns source_menu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from source_menu2


% --- Executes during object creation, after setting all properties.
function source_menu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_menu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in source_menu3.
function source_menu3_Callback(hObject, eventdata, handles)
kgrid = evalin('base','kgrid');
switch hObject.String{hObject.Value}
    case 'none'
        set(handles.source_text9,'Visible','off')
        set(handles.source_text10,'Visible','off')
        set(handles.source_text11,'Visible','off')
        set(handles.source_text12,'Visible','off')
        set(handles.source_mag3,'Visible','off')
        set(handles.source_xpos3,'Visible','off')
        set(handles.source_ypos3,'Visible','off')
        set(handles.source_rad3,'Visible','off')
         set(handles.sourcegroup3,'Visible','off');
    case 'Disk'
          set(handles.source_text9,'Visible','on')
        set(handles.source_text10,'Visible','on')
        set(handles.source_text11,'Visible','on')
        set(handles.source_text12,'Visible','on')
        set(handles.source_mag3,'Visible','on')
        set(handles.source_xpos3,'Visible','on')
        set(handles.source_ypos3,'Visible','on')
        set(handles.source_rad3,'Visible','on')
         set(handles.sourcegroup3,'Visible','on');
           if kgrid.Nz > 0
            set(handles.source_zpos3,'Visible','on')
        else
            set(handles.source_zpos3,'Visible','off')
        end
        
end
% hObject    handle to source_menu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns source_menu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from source_menu3


% --- Executes during object creation, after setting all properties.
function source_menu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_menu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_mag1_Callback(hObject, eventdata, handles)
% hObject    handle to source_mag1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_mag1 as text
%        str2double(get(hObject,'String')) returns contents of source_mag1 as a double


% --- Executes during object creation, after setting all properties.
function source_mag1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_mag1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_xpos1_Callback(hObject, eventdata, handles)
% hObject    handle to source_xpos1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_xpos1 as text
%        str2double(get(hObject,'String')) returns contents of source_xpos1 as a double


% --- Executes during object creation, after setting all properties.
function source_xpos1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_xpos1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_ypos1_Callback(hObject, eventdata, handles)
% hObject    handle to source_ypos1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_ypos1 as text
%        str2double(get(hObject,'String')) returns contents of source_ypos1 as a double


% --- Executes during object creation, after setting all properties.
function source_ypos1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_ypos1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_rad1_Callback(hObject, eventdata, handles)
% hObject    handle to source_rad1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_rad1 as text
%        str2double(get(hObject,'String')) returns contents of source_rad1 as a double


% --- Executes during object creation, after setting all properties.
function source_rad1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_rad1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_mag2_Callback(hObject, eventdata, handles)
% hObject    handle to source_mag2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_mag2 as text
%        str2double(get(hObject,'String')) returns contents of source_mag2 as a double


% --- Executes during object creation, after setting all properties.
function source_mag2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_mag2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_xpos2_Callback(hObject, eventdata, handles)
% hObject    handle to source_xpos2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_xpos2 as text
%        str2double(get(hObject,'String')) returns contents of source_xpos2 as a double


% --- Executes during object creation, after setting all properties.
function source_xpos2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_xpos2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_ypos2_Callback(hObject, eventdata, handles)
% hObject    handle to source_ypos2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_ypos2 as text
%        str2double(get(hObject,'String')) returns contents of source_ypos2 as a double


% --- Executes during object creation, after setting all properties.
function source_ypos2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_ypos2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_rad2_Callback(hObject, eventdata, handles)
% hObject    handle to source_rad2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_rad2 as text
%        str2double(get(hObject,'String')) returns contents of source_rad2 as a double


% --- Executes during object creation, after setting all properties.
function source_rad2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_rad2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_mag3_Callback(hObject, eventdata, handles)
% hObject    handle to source_mag3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_mag3 as text
%        str2double(get(hObject,'String')) returns contents of source_mag3 as a double


% --- Executes during object creation, after setting all properties.
function source_mag3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_mag3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_xpos3_Callback(hObject, eventdata, handles)
% hObject    handle to source_xpos3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_xpos3 as text
%        str2double(get(hObject,'String')) returns contents of source_xpos3 as a double


% --- Executes during object creation, after setting all properties.
function source_xpos3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_xpos3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_ypos3_Callback(hObject, eventdata, handles)
% hObject    handle to source_ypos3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_ypos3 as text
%        str2double(get(hObject,'String')) returns contents of source_ypos3 as a double


% --- Executes during object creation, after setting all properties.
function source_ypos3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_ypos3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_rad3_Callback(hObject, eventdata, handles)
% hObject    handle to source_rad3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_rad3 as text
%        str2double(get(hObject,'String')) returns contents of source_rad3 as a double


% --- Executes during object creation, after setting all properties.
function source_rad3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_rad3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function sensor_rad_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensor_rad_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in record_movie.
function record_movie_Callback(hObject, eventdata, handles)
% hObject    handle to record_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of record_movie



function sim_scale_Callback(hObject, eventdata, handles)
% hObject    handle to sim_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sim_scale as text
%        str2double(get(hObject,'String')) returns contents of sim_scale as a double


% --- Executes during object creation, after setting all properties.
function sim_scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sim_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sim_freq_Callback(hObject, eventdata, handles)
% hObject    handle to sim_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sim_freq as text
%        str2double(get(hObject,'String')) returns contents of sim_freq as a double


% --- Executes during object creation, after setting all properties.
function sim_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sim_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sim_fps_Callback(hObject, eventdata, handles)
% hObject    handle to sim_fps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sim_fps as text
%        str2double(get(hObject,'String')) returns contents of sim_fps as a double


% --- Executes during object creation, after setting all properties.
function sim_fps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sim_fps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sim_mesh.
function sim_mesh_Callback(hObject, eventdata, handles)
% hObject    handle to sim_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hit: get(hObject,'Value') returns toggle state of sim_mesh


% --- Executes on button press in sim_mask.
function sim_mask_Callback(hObject, eventdata, handles)
% hObject    handle to sim_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sim_mask


% --- Executes on button press in sim_pml.
function sim_pml_Callback(hObject, eventdata, handles)
% hObject    handle to sim_pml (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sim_pml


% --- Executes on button press in sensor_velocity.
function sensor_velocity_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_velocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sensor_velocity


% --- Executes on selection change in medium_menu.
function medium_menu_Callback(hObject, eventdata, handles)
V = hObject.Value;
switch V
    case 1
        set(handles.medium_layers,'Visible','off')
        set(handles.get_ct,'Visible','on')
    case 2
        set(handles.medium_layers,'Visible','on');
          set(handles.get_ct,'Visible','off')
end
% hObject    handle to medium_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns medium_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from medium_menu


% --- Executes during object creation, after setting all properties.
function medium_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to medium_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function medium_ss_Callback(hObject, eventdata, handles)
% hObject    handle to medium_ss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of medium_ss as text
%        str2double(get(hObject,'String')) returns contents of medium_ss as a double


% --- Executes during object creation, after setting all properties.
function medium_ss_CreateFcn(hObject, eventdata, handles)
% hObject    handle to medium_ss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function medium_density_Callback(hObject, eventdata, handles)
% hObject    handle to medium_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of medium_density as text
%        str2double(get(hObject,'String')) returns contents of medium_density as a double


% --- Executes during object creation, after setting all properties.
function medium_density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to medium_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function medium_alpha_coeff_Callback(hObject, eventdata, handles)
% hObject    handle to medium_alpha_coeff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of medium_alpha_coeff as text
%        str2double(get(hObject,'String')) returns contents of medium_alpha_coeff as a double


% --- Executes during object creation, after setting all properties.
function medium_alpha_coeff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to medium_alpha_coeff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function medium_alpha_power_Callback(hObject, eventdata, handles)
% hObject    handle to medium_alpha_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of medium_alpha_power as text
%        str2double(get(hObject,'String')) returns contents of medium_alpha_power as a double


% --- Executes during object creation, after setting all properties.
function medium_alpha_power_CreateFcn(hObject, eventdata, handles)
% hObject    handle to medium_alpha_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in data_menu.
function data_menu_Callback(hObject, eventdata, handles)
% hObject    handle to data_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns data_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from data_menu


% --- Executes during object creation, after setting all properties.
function data_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ly1_rnd.
function ly1_rnd_Callback(hObject, eventdata, handles)
% hObject    handle to ly1_rnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ly1_rnd



function ly3_yN_Callback(hObject, eventdata, handles)
% hObject    handle to ly3_yN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly3_yN as text
%        str2double(get(hObject,'String')) returns contents of ly3_yN as a double


% --- Executes during object creation, after setting all properties.
function ly3_yN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly3_yN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly2_yN_Callback(hObject, eventdata, handles)
% hObject    handle to ly2_yN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly2_yN as text
%        str2double(get(hObject,'String')) returns contents of ly2_yN as a double


% --- Executes during object creation, after setting all properties.
function ly2_yN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly2_yN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly1_yN_Callback(hObject, eventdata, handles)
% hObject    handle to ly1_yN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly1_yN as text
%        str2double(get(hObject,'String')) returns contents of ly1_yN as a double


% --- Executes during object creation, after setting all properties.
function ly1_yN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly1_yN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly3_y0_Callback(hObject, eventdata, handles)
% hObject    handle to ly3_y0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly3_y0 as text
%        str2double(get(hObject,'String')) returns contents of ly3_y0 as a double


% --- Executes during object creation, after setting all properties.
function ly3_y0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly3_y0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly3_xN_Callback(hObject, eventdata, handles)
% hObject    handle to ly3_xN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly3_xN as text
%        str2double(get(hObject,'String')) returns contents of ly3_xN as a double


% --- Executes during object creation, after setting all properties.
function ly3_xN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly3_xN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly3_x0_Callback(hObject, eventdata, handles)
% hObject    handle to ly3_x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly3_x0 as text
%        str2double(get(hObject,'String')) returns contents of ly3_x0 as a double


% --- Executes during object creation, after setting all properties.
function ly3_x0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly3_x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly2_y0_Callback(hObject, eventdata, handles)
% hObject    handle to ly2_y0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly2_y0 as text
%        str2double(get(hObject,'String')) returns contents of ly2_y0 as a double


% --- Executes during object creation, after setting all properties.
function ly2_y0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly2_y0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly2_xN_Callback(hObject, eventdata, handles)
% hObject    handle to ly2_xN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly2_xN as text
%        str2double(get(hObject,'String')) returns contents of ly2_xN as a double


% --- Executes during object creation, after setting all properties.
function ly2_xN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly2_xN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly2_x0_Callback(hObject, eventdata, handles)
% hObject    handle to ly2_x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly2_x0 as text
%        str2double(get(hObject,'String')) returns contents of ly2_x0 as a double


% --- Executes during object creation, after setting all properties.
function ly2_x0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly2_x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly1_y0_Callback(hObject, eventdata, handles)
% hObject    handle to ly1_y0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly1_y0 as text
%        str2double(get(hObject,'String')) returns contents of ly1_y0 as a double


% --- Executes during object creation, after setting all properties.
function ly1_y0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly1_y0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly1_xN_Callback(hObject, eventdata, handles)
% hObject    handle to ly1_xN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly1_xN as text
%        str2double(get(hObject,'String')) returns contents of ly1_xN as a double


% --- Executes during object creation, after setting all properties.
function ly1_xN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly1_xN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly1_x0_Callback(hObject, eventdata, handles)
% hObject    handle to ly1_x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly1_x0 as text
%        str2double(get(hObject,'String')) returns contents of ly1_x0 as a double


% --- Executes during object creation, after setting all properties.
function ly1_x0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly1_x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ly1_on.
function ly1_on_Callback(hObject, eventdata, handles)
% hObject    handle to ly1_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ly1_on


% --- Executes on button press in ly2_on.
function ly2_on_Callback(hObject, eventdata, handles)
% hObject    handle to ly2_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ly2_on


% --- Executes on button press in ly3_on.
function ly3_on_Callback(hObject, eventdata, handles)
% hObject    handle to ly3_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ly3_on



function edit47_Callback(hObject, eventdata, handles)
% hObject    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit47 as text
%        str2double(get(hObject,'String')) returns contents of edit47 as a double


% --- Executes during object creation, after setting all properties.
function edit47_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly3_ac_Callback(hObject, eventdata, handles)
% hObject    handle to ly3_ac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly3_ac as text
%        str2double(get(hObject,'String')) returns contents of ly3_ac as a double


% --- Executes during object creation, after setting all properties.
function ly3_ac_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly3_ac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly3_den_Callback(hObject, eventdata, handles)
% hObject    handle to ly3_den (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly3_den as text
%        str2double(get(hObject,'String')) returns contents of ly3_den as a double


% --- Executes during object creation, after setting all properties.
function ly3_den_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly3_den (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly3_ss_Callback(hObject, eventdata, handles)
% hObject    handle to ly3_ss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly3_ss as text
%        str2double(get(hObject,'String')) returns contents of ly3_ss as a double


% --- Executes during object creation, after setting all properties.
function ly3_ss_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly3_ss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit43_Callback(hObject, eventdata, handles)
% hObject    handle to edit43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit43 as text
%        str2double(get(hObject,'String')) returns contents of edit43 as a double


% --- Executes during object creation, after setting all properties.
function edit43_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly2_ac_Callback(hObject, eventdata, handles)
% hObject    handle to ly2_ac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly2_ac as text
%        str2double(get(hObject,'String')) returns contents of ly2_ac as a double


% --- Executes during object creation, after setting all properties.
function ly2_ac_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly2_ac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly2_den_Callback(hObject, eventdata, handles)
% hObject    handle to ly2_den (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly2_den as text
%        str2double(get(hObject,'String')) returns contents of ly2_den as a double


% --- Executes during object creation, after setting all properties.
function ly2_den_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly2_den (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly2_ss_Callback(hObject, eventdata, handles)
% hObject    handle to ly2_ss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly2_ss as text
%        str2double(get(hObject,'String')) returns contents of ly2_ss as a double


% --- Executes during object creation, after setting all properties.
function ly2_ss_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly2_ss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit39_Callback(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit39 as text
%        str2double(get(hObject,'String')) returns contents of edit39 as a double


% --- Executes during object creation, after setting all properties.
function edit39_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly1_ac_Callback(hObject, eventdata, handles)
% hObject    handle to ly1_ac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly1_ac as text
%        str2double(get(hObject,'String')) returns contents of ly1_ac as a double


% --- Executes during object creation, after setting all properties.
function ly1_ac_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly1_ac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly1_den_Callback(hObject, eventdata, handles)
% hObject    handle to ly1_den (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly1_den as text
%        str2double(get(hObject,'String')) returns contents of ly1_den as a double


% --- Executes during object creation, after setting all properties.
function ly1_den_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly1_den (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly1_ss_Callback(hObject, eventdata, handles)
% hObject    handle to ly1_ss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly1_ss as text
%        str2double(get(hObject,'String')) returns contents of ly1_ss as a double


% --- Executes during object creation, after setting all properties.
function ly1_ss_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly1_ss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_medium.
function plot_medium_Callback(hObject, eventdata, handles)
medium = evalin('base','medium');
medium.sound_speed = flipud(medium.sound_speed);
medium.alpha_coeff = flipud(medium.alpha_coeff);
medium.density = flipud(medium.density);
kgrid = evalin('base','kgrid');
if length(medium.sound_speed) == 1
    medium.sound_speed = ones(kgrid.Ny,kgrid.Nx)*medium.sound_speed;
end
x = (kgrid.y_vec-min(kgrid.y_vec))*100;
y = (kgrid.x_vec-min(kgrid.y_vec))*100;
y = flipud(y);
if handles.ef.Value
    figure(3);
else
axes(handles.axes2)
end
imagesc(x,y,medium.sound_speed);
xlabel('Depth (cm)');
ylabel('Lateral (cm)');
title('Sound Speed');
% colormap('jet');
colormap(ones(64,3));
colormap('gray');
colorbar;

if handles.ef.Value
    figure(4);
else
axes(handles.axes3)
end
hold off;
if length(medium.density) == 1
    medium.density = ones(kgrid.Ny,kgrid.Nx)*medium.density;
end
imagesc(x,y,medium.density);
xlabel('Depth (cm)');
ylabel('Lateral (cm)');
title('Density');
colormap('gray');
colorbar;

if handles.ef.Value
    figure(5);
else
axes(handles.axes1)
end
hold off;
if length(medium.alpha_coeff) == 1
    medium.density = ones(kgrid.Ny,kgrid.Nx)*medium.alpha_coeff;
end
imagesc(x,y,medium.alpha_coeff);
xlabel('Depth (cm)');
ylabel('Lateral (cm)');
title('Alpha Coeff');
colormap('gray');
colorbar;


% hObject    handle to plot_medium (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function source_delay1_Callback(hObject, eventdata, handles)
% hObject    handle to source_delay1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_delay1 as text
%        str2double(get(hObject,'String')) returns contents of source_delay1 as a double


% --- Executes during object creation, after setting all properties.
function source_delay1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_delay1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_duration1_Callback(hObject, eventdata, handles)
% hObject    handle to source_duration1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_duration1 as text
%        str2double(get(hObject,'String')) returns contents of source_duration1 as a double


% --- Executes during object creation, after setting all properties.
function source_duration1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_duration1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in source_type1.
function source_type1_Callback(hObject, eventdata, handles)
% hObject    handle to source_type1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns source_type1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from source_type1


% --- Executes during object creation, after setting all properties.
function source_type1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_type1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_period1_Callback(hObject, eventdata, handles)
% hObject    handle to source_period1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_period1 as text
%        str2double(get(hObject,'String')) returns contents of source_period1 as a double


% --- Executes during object creation, after setting all properties.
function source_period1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_period1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sensor_final.
function sensor_final_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_final (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sensor_final



function display_clims_Callback(hObject, eventdata, handles)
% hObject    handle to display_clims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of display_clims as text
%        str2double(get(hObject,'String')) returns contents of display_clims as a double


% --- Executes during object creation, after setting all properties.
function display_clims_CreateFcn(hObject, eventdata, handles)
% hObject    handle to display_clims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in trans_menu.
function trans_menu_Callback(hObject, eventdata, handles)
trans = handles.trans_menu.String{handles.trans_menu.Value};
kgrid = evalin('base','kgrid');
switch trans
    case '1MHz'
        set(handles.source_period1,'String',1)
        set(handles.source_duration1,'String',2)
        set(handles.source_delay1,'String',0)
        x = kgrid.Nx/2;
        y = kgrid.Ny-4;
        set(handles.trans_arc_pos,'String',num2str([x y]))
        r = 68/kgrid.dx/1e3;
        d = round((30/kgrid.dx)/1e3)+1;
        set(handles.trans_radius,'String',r);
        set(handles.trans_diameter,'String',d);
        f = [x round(68/kgrid.dx/1e3)-4];
        set(handles.trans_focus,'String',num2str(f));
        set(handles.arraygroup,'Visible','off');
        set(handles.transwave,'Visible','off');
           set(handles.trans_duration,'Visible','on');
        set(handles.trans_period,'Visible','on');
        set(handles.trans_type,'Visible','on');
        set(handles.trans_arc_pos,'Visible','on');
        set(handles.trans_radius,'Visible','on');
        set(handles.trans_diameter,'Visible','on');
    case 'P4-1'
        x = kgrid.Nx/2;
        set(handles.trans_x0,'String',x); %grid
        y = 10;
        set(handles.trans_y0,'String',y); %grid
        Trans.kerf = .051;
        Trans.width = .2443;
        Trans.space = Trans.kerf+Trans.width;
        %         xkerf = 1/kgrid.dx/1e3;
        %         length = 1/kgrid.dx/1e3;
        dx = kgrid.dx*1000;
        kerf = round(Trans.kerf/dx);
        width = round(Trans.width/dx);
        set(handles.trans_x_kerf,'String',num2str(kerf)); %mm
        set(handles.trans_width,'String',1);
        set(handles.trans_length,'String',width);
        set(handles.trans_cycles,'String',3);
        set(handles.trans_freq,'String',2.5); %MHz
        set(handles.trans_steer,'String',0); %degrees
        set(handles.trans_el_focus,'String',80); %mm
        set(handles.trans_lens_focus,'String',60) %mm
        set(handles.trans_elements,'String',96)
        set(handles.trans_active_elements,'String','1:96');
        set(handles.arraygroup,'Visible','on');
        set(handles.transwave,'Visible','on');
        set(handles.trans_duration,'Visible','off');
        set(handles.trans_period,'Visible','off');
        set(handles.trans_type,'Visible','off');
        set(handles.trans_arc_pos,'Visible','off');
        set(handles.trans_radius,'Visible','off');
        set(handles.trans_diameter,'Visible','off');
         set(handles.trans_focus,'Visible','on');
    case 'P4-2'
        x = kgrid.Nx/2;
        set(handles.trans_x0,'String',x); %grid
        y = 10;
        set(handles.trans_y0,'String',y); %grid
        xkerf = 0.3/kgrid.dx/1e3;
        length = 0.4831/kgrid.dx/1e3;
        set(handles.trans_x_kerf,'String',num2str(xkerf)); %mm
        set(handles.trans_width,'String',num2str(length));
        set(handles.trans_length,'String',num2str(length));
        set(handles.trans_cycles,'String',3);
        set(handles.trans_freq,'String',3); %MHz
        set(handles.trans_steer,'String',0); %degrees
        set(handles.trans_el_focus,'String',60); %mm
        set(handles.trans_lens_focus,'String',60) %mm
        set(handles.trans_active_elements,'String','1:64');
        set(handles.arraygroup,'Visible','on');
        set(handles.transwave,'Visible','on');
    case 'H235'
        set(handles.arraygroup,'Visible','on');
        set(handles.transwave,'Visible','on');
    case 'h247'
        set(handles.arraygroup,'Visible','on');
        set(handles.transwave,'Visible','on');
    case 'Custom'
        x = kgrid.Nx/2;
        set(handles.trans_x0,'String',x); %grid
        y = 10;
        set(handles.trans_y0,'String',y); %grid
        xkerf = 1/kgrid.dx/1e3-1;
        length = 1;
        set(handles.trans_x_kerf,'String',num2str(xkerf)); %mm
        set(handles.trans_width,'String',num2str(1));
        set(handles.trans_length,'String',num2str(length));
        set(handles.trans_cycles,'String',3);
        set(handles.trans_freq,'String',1); %MHz
        set(handles.trans_steer,'String',0); %degrees
        set(handles.trans_el_focus,'String',60); %mm
        set(handles.trans_lens_focus,'String',60) %mm
        set(handles.trans_elements,'String',25)
        set(handles.trans_active_elements,'String','1:25');
        set(handles.arraygroup,'Visible','on');
        set(handles.transwave,'Visible','on');
        set(handles.trans_duration,'Visible','off');
        set(handles.trans_period,'Visible','off');
        set(handles.trans_type,'Visible','off');
        set(handles.trans_arc_pos,'Visible','off');
        set(handles.trans_radius,'Visible','off');
        set(handles.trans_diameter,'Visible','off');
end
% hObject    handle to trans_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns trans_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from trans_menu


% --- Executes during object creation, after setting all properties.
function trans_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_x0_Callback(hObject, eventdata, handles)
% hObject    handle to trans_x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_x0 as text
%        str2double(get(hObject,'String')) returns contents of trans_x0 as a double


% --- Executes during object creation, after setting all properties.
function trans_x0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_y0_Callback(hObject, eventdata, handles)
% hObject    handle to trans_y0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_y0 as text
%        str2double(get(hObject,'String')) returns contents of trans_y0 as a double


% --- Executes during object creation, after setting all properties.
function trans_y0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_y0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_arc_pos_Callback(hObject, eventdata, handles)
% hObject    handle to trans_arc_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_arc_pos as text
%        str2double(get(hObject,'String')) returns contents of trans_arc_pos as a double


% --- Executes during object creation, after setting all properties.
function trans_arc_pos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_arc_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_radius_Callback(hObject, eventdata, handles)
% hObject    handle to trans_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_radius as text
%        str2double(get(hObject,'String')) returns contents of trans_radius as a double


% --- Executes during object creation, after setting all properties.
function trans_radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_diameter_Callback(hObject, eventdata, handles)
% hObject    handle to trans_diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_diameter as text
%        str2double(get(hObject,'String')) returns contents of trans_diameter as a double


% --- Executes during object creation, after setting all properties.
function trans_diameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_focus_Callback(hObject, eventdata, handles)
% hObject    handle to trans_focus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_focus as text
%        str2double(get(hObject,'String')) returns contents of trans_focus as a double


% --- Executes during object creation, after setting all properties.
function trans_focus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_focus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function record_filename_Callback(hObject, eventdata, handles)
% hObject    handle to record_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of record_filename as text
%        str2double(get(hObject,'String')) returns contents of record_filename as a double


% --- Executes during object creation, after setting all properties.
function record_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to record_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function record_dir_Callback(hObject, eventdata, handles)
% hObject    handle to record_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of record_dir as text
%        str2double(get(hObject,'String')) returns contents of record_dir as a double


% --- Executes during object creation, after setting all properties.
function record_dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to record_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_elements_Callback(hObject, eventdata, handles)
% hObject    handle to trans_elements (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_elements as text
%        str2double(get(hObject,'String')) returns contents of trans_elements as a double


% --- Executes during object creation, after setting all properties.
function trans_elements_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_elements (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_x_kerf_Callback(hObject, eventdata, handles)
% hObject    handle to trans_x_kerf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_x_kerf as text
%        str2double(get(hObject,'String')) returns contents of trans_x_kerf as a double


% --- Executes during object creation, after setting all properties.
function trans_x_kerf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_x_kerf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_steer_Callback(hObject, eventdata, handles)
% hObject    handle to trans_steer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_steer as text
%        str2double(get(hObject,'String')) returns contents of trans_steer as a double


% --- Executes during object creation, after setting all properties.
function trans_steer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_steer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_freq_Callback(hObject, eventdata, handles)
% hObject    handle to trans_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_freq as text
%        str2double(get(hObject,'String')) returns contents of trans_freq as a double


% --- Executes during object creation, after setting all properties.
function trans_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_cycles_Callback(hObject, eventdata, handles)
% hObject    handle to trans_cycles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_cycles as text
%        str2double(get(hObject,'String')) returns contents of trans_cycles as a double


% --- Executes during object creation, after setting all properties.
function trans_cycles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_cycles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sim_3D.
function sim_3D_Callback(hObject, eventdata, handles)
% hObject    handle to sim_3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sim_3D


% --- Executes on button press in make_trans.
function make_trans_Callback(hObject, eventdata, handles)
kgrid = evalin('base','kgrid');
medium = evalin('base','medium');
% tw = evalin('base','tw');
transducer.number_elements = str2double(handles.trans_elements.String);    % total number of transducer elements
e_width = str2double(handles.trans_width.String); 
e_width = round(e_width/kgrid.dy/1e3);
transducer.element_width = e_width;       % width of each element [grid points]
e_length = str2double(handles.trans_length.String);
e_length = round(e_length/kgrid.dz/1e3);
transducer.element_length = e_length;    % length of each element [grid points]
e_kerf = str2double(handles.trans_x_kerf.String);
e_kerf = round(e_kerf/kgrid.dy/1e3);
transducer.element_spacing = e_kerf;    % spacing (kerf width) between the elements [grid points]
transducer.radius = inf;           % radius of curvature of the transducer [m]
% Denotes center of transducer position
x0 = str2double(handles.trans_x0.String);
y0 = str2double(handles.trans_y0.String);
z0 = str2double(handles.trans_z0.String);
%Do I need to normalize transducer width to grid elements still?
transducer.position = round([x0,y0,z0]);

%Beamforming delay properities
%transducer.sound_speed = medium.sound_speed; %[m/s]
transducer.sound_speed = medium.sound_speed(1,1); %[m/s]
% focus_loc = str2num(handles.trans_focus.String);
% focus = [kgrid.x_vec(focus_loc(1)) kgrid.y_vec(focus_loc(2)) kgrid.z_vec(focus_loc(3))];
transducer.focus_distance = str2double(handles.trans_lens_focus.String)/1e3; %[m]
transducer.elevation_focus_distance = str2double(handles.trans_el_focus.String)/1e3; % [m]
transducer.steering_angle = str2double(handles.trans_steer.String); %[deg]

% Apodization and active elements
transducer.transmit_apodization = handles.t_apod.String{handles.t_apod.Value};
transducer.receive_apodization = handles.r_apod.String{handles.r_apod.Value};
transducer.active_elements = zeros(transducer.number_elements,1);
active = str2num(handles.trans_active_elements.String);
transducer.active_elements(active) = 1;

% Create transducer
tw = evalin('base','tw');
transducer.input_signal = tw;
transducer = kWaveTransducer(kgrid, transducer);
assignin('base','transducer',transducer);


% hObject    handle to make_trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in show_trans.
function show_trans_Callback(hObject, eventdata, handles)
transducer = evalin('base','transducer');
voxelPlot(double(transducer.all_elements_mask));
% hObject    handle to show_trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in show_tw.
function show_tw_Callback(hObject, eventdata, handles)
kgrid = evalin('base','kgrid');
medium = evalin('base','medium');
strength = str2double(handles.trans_mag.String)*1e6;
freq = str2double(handles.trans_freq.String)*1e6;
cycles = str2double(handles.trans_cycles.String);
envel = handles.tw_encoding.String{handles.tw_encoding.Value};
switch envel
    case 'Gauss'
        signal = toneBurst(1/kgrid.dt,freq,cycles);
    case 'Rect'
        signal = toneBurst(1/kgrid.dt,freq,cycles,'Envelope','Rectangular');
end
signal2 = interp1(linspace(0,1,length(signal)),signal,linspace(0,1,length(signal)*10));
signal = signal.*strength;
signal2 = signal2*strength;
t = interp1(linspace(0,1,length(kgrid.t_array)),kgrid.t_array,linspace(0,1,length(kgrid.t_array)*10));
if handles.ef.Value
    figure(555)
%     signal = resample(signal,10,1);
%     t_array = resample(kgrid.t_array,10,1);
%     t_array = t_array - min(t_array);
    plot(t(1:length(signal)),signal,'k','LineWidth',2.5);
    xlabel('Time [us]')
    ylabel('Particle Velocity [m/s]')
    
    figure(556)
    dt = kgrid.dt;
    maxf = 1/dt/2;
%     signal = interp1(1:length(signal),signal,1:0.25:length(signal));
    f_axis = linspace(0,maxf,length(signal2)*4);
    Signal = fft(padarray(signal2',length(signal2)*4-round(length(signal2)/2),'both'));
%     Signal = interp1(1:length(Signal),Signal,1:length(Signal)*10);
    plot(f_axis/1e6,abs(Signal(1:length(f_axis))),'k','LineWidth',2.5)
    xlim([0 10])
    xlabel('Frequency [MHz]')
    ylabel('Amplitude [au]')    
else
    axes(handles.axes2)
    %plot(t(1:length(signal2)),signal2,'k','LineWidth',2.5);
    plot(kgrid.t_array(1:length(signal))*1e6,signal);
    xlabel('Time [us]')
    ylabel('Particle Velocity [m/s]')
    axes(handles.axes3)
    dt = kgrid.dt;
    maxf = 1/dt/2;
   % f_axis = linspace(0,maxf,length(signal2)/2);
     f_axis = linspace(0,maxf,length(signal)/2);
    Signal = fft(signal);
    Signal2 = interp1(linspace(0,1,length(signal)),abs(Signal),linspace(0,1,length(signal2)));
    plot(f_axis/1e6,abs(Signal(1:length(f_axis))),'k','LineWidth',2.5)
  %  plot(f_axis/1e6,abs(Signal2(1:length(f_axis))),'k','LineWidth',2.5)
    xlim([0 10])
    xlabel('Frequency [MHz]')
    ylabel('Amplitude [au]')
end
assignin('base','tw',signal2);

% hObject    handle to show_tw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function trans_width_Callback(hObject, eventdata, handles)
% hObject    handle to trans_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_width as text
%        str2double(get(hObject,'String')) returns contents of trans_width as a double


% --- Executes during object creation, after setting all properties.
function trans_width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_length_Callback(hObject, eventdata, handles)
% hObject    handle to trans_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_length as text
%        str2double(get(hObject,'String')) returns contents of trans_length as a double


% --- Executes during object creation, after setting all properties.
function trans_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in wave_type.
function wave_type_Callback(hObject, eventdata, handles)
% if hObject.Value == 1
%      set(handles.hadamard,'String','Planes')
%      set(handles.text138,'String','dTheta (deg)');
% elseif hObject.Value == 3
%     set(handles.hadamard,'String','Hadamard')
%      set(handles.text138,'String','STA');
% end
% hObject    handle to wave_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns wave_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from wave_type


% --- Executes during object creation, after setting all properties.
function wave_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wave_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in center_trans.
function center_trans_Callback(hObject, eventdata, handles)
kgrid = evalin('base','kgrid')
y = kgrid.Ny/2;
x = kgrid.Nx/2;
set(handles.trans_x0,'String',x)
set(handles.trans_y0,'String',y)
if kgrid.Nz > 0
    z = kgrid.Nz/2;
    set(handles.trans_z0,'String',z);
end
% hObject    handle to center_trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function trans_z0_Callback(hObject, eventdata, handles)
% hObject    handle to trans_z0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_z0 as text
%        str2double(get(hObject,'String')) returns contents of trans_z0 as a double


% --- Executes during object creation, after setting all properties.
function trans_z0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_z0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_lens_focus_Callback(hObject, eventdata, handles)
% hObject    handle to trans_lens_focus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_lens_focus as text
%        str2double(get(hObject,'String')) returns contents of trans_lens_focus as a double


% --- Executes during object creation, after setting all properties.
function trans_lens_focus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_lens_focus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_el_focus_Callback(hObject, eventdata, handles)
% hObject    handle to trans_el_focus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_el_focus as text
%        str2double(get(hObject,'String')) returns contents of trans_el_focus as a double


% --- Executes during object creation, after setting all properties.
function trans_el_focus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_el_focus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in t_apod.
function t_apod_Callback(hObject, eventdata, handles)
% hObject    handle to t_apod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns t_apod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from t_apod


% --- Executes during object creation, after setting all properties.
function t_apod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_apod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in r_apod.
function r_apod_Callback(hObject, eventdata, handles)
% hObject    handle to r_apod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns r_apod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from r_apod


% --- Executes during object creation, after setting all properties.
function r_apod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r_apod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_active_elements_Callback(hObject, eventdata, handles)
% hObject    handle to trans_active_elements (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_active_elements as text
%        str2double(get(hObject,'String')) returns contents of trans_active_elements as a double


% --- Executes during object creation, after setting all properties.
function trans_active_elements_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_active_elements (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sensor_zpos_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_zpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensor_zpos as text
%        str2double(get(hObject,'String')) returns contents of sensor_zpos as a double


% --- Executes during object creation, after setting all properties.
function sensor_zpos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensor_zpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_mag_Callback(hObject, eventdata, handles)
% hObject    handle to trans_mag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_mag as text
%        str2double(get(hObject,'String')) returns contents of trans_mag as a double


% --- Executes during object creation, after setting all properties.
function trans_mag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_mag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_delay_Callback(hObject, eventdata, handles)
% hObject    handle to trans_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_delay as text
%        str2double(get(hObject,'String')) returns contents of trans_delay as a double


% --- Executes during object creation, after setting all properties.
function trans_delay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_duration_Callback(hObject, eventdata, handles)
% hObject    handle to trans_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_duration as text
%        str2double(get(hObject,'String')) returns contents of trans_duration as a double


% --- Executes during object creation, after setting all properties.
function trans_duration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_period_Callback(hObject, eventdata, handles)
% hObject    handle to trans_period (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_period as text
%        str2double(get(hObject,'String')) returns contents of trans_period as a double


% --- Executes during object creation, after setting all properties.
function trans_period_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_period (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in trans_type.
function trans_type_Callback(hObject, eventdata, handles)
% hObject    handle to trans_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns trans_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from trans_type


% --- Executes during object creation, after setting all properties.
function trans_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in tw_encoding.
function tw_encoding_Callback(hObject, eventdata, handles)
% hObject    handle to tw_encoding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tw_encoding contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tw_encoding


% --- Executes during object creation, after setting all properties.
function tw_encoding_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tw_encoding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_zpos1_Callback(hObject, eventdata, handles)
% hObject    handle to source_zpos1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_zpos1 as text
%        str2double(get(hObject,'String')) returns contents of source_zpos1 as a double


% --- Executes during object creation, after setting all properties.
function source_zpos1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_zpos1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_zpos3_Callback(hObject, eventdata, handles)
% hObject    handle to source_zpos3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_zpos3 as text
%        str2double(get(hObject,'String')) returns contents of source_zpos3 as a double


% --- Executes during object creation, after setting all properties.
function source_zpos3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_zpos3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function source_zpos2_Callback(hObject, eventdata, handles)
% hObject    handle to source_zpos2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_zpos2 as text
%        str2double(get(hObject,'String')) returns contents of source_zpos2 as a double


% --- Executes during object creation, after setting all properties.
function source_zpos2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_zpos2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sensor_p_max.
function sensor_p_max_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_p_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sensor_p_max


% --- Executes on button press in sensor_p_max_all.
function sensor_p_max_all_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_p_max_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sensor_p_max_all


% --- Executes on button press in custom_sensor.
function custom_sensor_Callback(hObject, eventdata, handles)
ksim_sensor;
% hObject    handle to custom_sensor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in get_ct.
function get_ct_Callback(hObject, eventdata, handles)
load_ct;
% hObject    handle to get_ct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
% bg = ones(str2double(handles.medx.String), str2double(handles.medy.String));
% med.k = bg*str2double(handles.medk.String);
% med.rho = bg*str2double(handles.medrho.String);
% x = str2num(handles.elecx.String);
% y = str2num(handles.elecy.String);


% Get Pressure
psource = handles.pressuremenu.String{handles.pressuremenu.Value};
switch psource
%     case 'sensor_data'
%         sensor_data = evalin('base','sensor_data');
%         pressure = sensor_data.p;
%         kgrid = evalin('base','kgrid');
%         pressure = reshape(pressure,[kgrid.Nx,kgrid.Ny,size(pressure,2)]);
%         assignin('base','pressure',pressure);
    case 'pressure'
        pressure = evalin('base','pressure');
    case 'convolution'
        P_all = evalin('base','P_all');
    case 'multi_conv'
        r = str2double(get(handles.auto_input1,'String'));
        load('D:\Chet Backup\Manuscripts\Modeling\Russp1','delays')
    case {'Hadamard', 'planes', 'H_pres'}
        P = evalin('base','H_pressure');
end


% Get Current
kgrid = evalin('base','kgrid');
sensor = evalin('base','sensor');
if ~handles.use_elec.Value %Use the ksim electrode design
    bottomleft = [find(kgrid.x_vec == min(sensor.mask(1,:)),1) find(kgrid.y_vec == min(sensor.mask(2,:)),1)];
    I_amp = str2double(handles.camp.String);
    I_freq = str2double(handles.cfreq.String);
    I_dur = str2double(handles.cdur.String);
    I_dt = str2double(handles.cdt.String);
    ctime = linspace(0,I_dur,I_dur/I_dt)/1000;
    I_delay = str2double(handles.cdelay.String);
    I = zeros(v,v);
    I_loc = [x - bottomleft(1), y - bottomleft(2)];
    med.k(I_loc(1),I_loc(2)) = str2double(handles.eleck.String);
    med.rho(I_loc(1),I_loc(2)) = str2double(handles.elecrho.String);

    if I_freq > 0
        I_freq = I_freq;
        I_mod = cos(2*pi*I_freq*ctime);
        I_amp2 = I_amp.*I_mod;
    else
        I_amp2(1:length(ctime)) = I_amp;
    end
    I(I_loc(1),I_loc(2),1:length(ctime)) = I_amp2;

    % NEED CURRENT SPREAD INTO SURROUNDING TISSUE


    % b = gausswin(7);
    %Simple method for approximating flow of current into surrounding tissue
    for i = 1:size(I,3)
        x2(:,:,i) = imgaussfilt(I(:,:,i),1.7);
        %             x2(i,j,k,:) = filter(b,1,squeeze(x(i,j,k,:)));
    end
    I = x2;
    clear x2;

    I = repmat(I,[1,1,1,kgrid.Nt]);
    I = permute(I,[1,2,4,3]);
    assignin('base','I',I);

else % Import from Electroid GUI
    electroids = evalin('base','electroids');
    med.k = electroids.k(1);
    med.rho = electroids.rho(1);
    I = evalin('base','I');
    I = I(:,:,1,1);
    ctime = 1;
    I = repmat(I,[1,1,kgrid.Nt]);
    %     I = permute(I,[1,2,4,3]);
end

%Simple method for amplifying surrounding tissue based on conductivity
% for i = 1:size(I,3)
%     I(:,:,i) = I(:,:,i).*med.rho;
% end

%Importing and using lead field
if handles.uself.Value

    L = evalin('base','L');
    %     for i = 1:size(I,3)
    %         I(:,:,i) = I(:,:,i).*L;
    %     end

end

% Deducing dot product from gradients




% Get V by multiplying pressure matrix by current matrix over time
% V = ones(size(I),'gpuArray');
switch psource
    case 'convolution'

        delays = evalin('base','delays');
        V = zeros(length(P_all),size(P_all{1},ndims(P_all{1})));
        for i = 1:length(P_all)
            for j = 1:size(P_all{1},ndims(P_all{1}))
                V(i,j) = sum(sum(P_all{i}(:,:,j).*I(:,:,j,1).*-med.k.*med.rho));
                %             V(i,j) = P_all{i}(:,:,j),I(:,:,j,1);
            end
            V(i,:) = circshift(V(i,:),-delays(i));
            multiWaitbar('Solving Forward Solution',i/length(P_all));
        end
    case {'sensor_data','pressure'}

        V = ones([size(I)]);
        I = I(:,:,1);
        if handles.uself.Value
            [Ix Iy] = gradient(I);
            for k = 1:size(L,3)
                [Lx(:,:,k) Ly(:,:,k)] = gradient(L(:,:,k));
            end
            for i = 1:size(L,1)
                for j = 1:size(L,2)
                    for k = 1:size(L,3)
                        LFM(i,j,k) = Ix(i,j)*Lx(i,j,k)+Iy(i,j)*Ly(i,j,k);
                    end
                end
                multiWaitbar('Generating Lead Field',i/size(L,1));
            end
        end
        for i = 1:kgrid.Nt
            for j = 1:length(ctime)
                if handles.uself.Value
                    for k = 1:size(L,3)
                        %                     V(:,:,i,j,k) = pressure(:,:,i).*(I(:,:,i,j).*L(:,:,k)).*med.k.*med.rho;
                        V(:,:,i,j,k) = pressure(:,:,i).*I(:,:,j).*LFM(:,:,k).*med.k.*med.rho;
                        M(i,k,:) = squeeze(sum(sum(V(:,:,i,j,k))));
                    end
                else
                    V(:,:,i,j) = pressure(:,:,i).*I(:,:,j).*med.k.*med.rho;
                    M(i) = squeeze(sum(sum(V(:,:,i,j))));
                end
            end
            multiWaitbar('Solving Forward Solution',i/kgrid.Nt);
        end
        assignin('base','M',M);
    case 'multi_conv'
        V = zeros(length(P_all)*2,size(P_all{1},3));
        for i = 1:length(P_all)
            for j = 1:size(P_all{1},3)
                V(:,:,i,j)= sum(conv2(P_all{i}(:,:,j),I(:,:,j,1).*-med.k.*med.rho));
            end %There should be a dividing factor somewhere
            V(i,:) = circshift(V(i,:),-delays(i));
            multiWaitbar('Solving Forward Solution',i/length(P_all));
        end
        clear P_all
        clear delays
        load('D:\Chet Backup\Manuscripts\Modeling\Russp2',['P_all','delays'])
        for i = 1:length(P_all)
            for j = 1:size(P_all{1},3)
                V(:,:,i+length(P_all),j)= sum(conv2(P_all{i}(:,:,j),I(:,:,j,1).*-med.k.*med.rho));
            end %There should be a dividing factor somewhere
            V(i+length(P_all),:) = circshift(V(i+length(P_all),:),-delays(i));
            multiWaitbar('Solving Forward Solution',i/length(P_all));
        end
    case {'Hadamard'} %%% This only exists for 2 dimensions for now
        V2 = ones(size(P));
        I = I(:,:,1);
        if handles.uself.Value
            [Ix Iy] = gradient(I);
            for k = 1:size(L,3)
                [Lx(:,:,k) Ly(:,:,k)] = gradient(L(:,:,k));
            end
            for i = 1:size(L,1)
                for j = 1:size(L,2)
                    for k = 1:size(L,3)
                        LFM(i,j,k) = Ix(i,j)*Lx(i,j,k)+Iy(i,j)*Ly(i,j,k);
                    end
                end
                multiWaitbar('Generating Lead Field',i/size(L,1));
            end
        end
        for i = 1:kgrid.Nt
            for j = 1:length(ctime)
                for m = 1:size(P,3)
                    if handles.uself.Value
                        for k = 1:size(L,3)
                            %                     V(:,:,i,j,k) = pressure(:,:,i).*(I(:,:,i,j).*L(:,:,k)).*med.k.*med.rho;
                            V2(:,:,i,j,k) = P(:,:,i).*I(:,:,j).*LFM(:,:,k).*med.k.*med.rho;
                            M(i,k,:) = squeeze(sum(sum(V(:,:,i,j,k))));
                        end
                    else
                        V2(:,:,m,i,j) = P(:,:,m,i).*I(:,:,j).*med.k.*med.rho;
                        M(i) = squeeze(sum(sum(V2(:,:,m,i,j))));
                    end
                end
            end
            multiWaitbar('Solving Forward Solution',i/kgrid.Nt);
        end
        V2 = squeeze(V2);

        xmax = max(V2(:)).*2; %This *2 is to make up for the -.5 in the next var
        %    noises = (rand([size(V2,1),size(V2,2)])-.5).*xmax; % Too large
        snrs = str2double(handles.noise_snr.String);
        %  noise_var = noises./snrs;
        for i = 1:size(V2,3)
            for j = 1:size(V2,4)
                V2(:,:,i,j) = V2(:,:,i,j) +  (rand([size(V2,1),size(V2,2)])-.5).*xmax./snrs;
            end
            multiWaitbar('Adding Noise',i/size(V2,3));
        end
        %Xnoise = V2 + noise_var;

        V = squeeze(sum(V2,3));
        multiWaitbar('CLOSEALL');
    case 'H_pres'
        V2 = ones(size(P));
        I = I(:,:,1);
        if handles.uself.Value
            [Ix Iy] = gradient(I);
            for k = 1:size(L,3)
                [Lx(:,:,k) Ly(:,:,k)] = gradient(L(:,:,k));
            end
            for i = 1:size(L,1)
                for j = 1:size(L,2)
                    for k = 1:size(L,3)
                        LFM(i,j,k) = Ix(i,j)*Lx(i,j,k)+Iy(i,j)*Ly(i,j,k);
                    end
                end
                multiWaitbar('Generating Lead Field',i/size(L,1));
            end
        end

        %%%% Add Noise
        I2 = I;
        for m = 1:size(P,3) %Transmit number
            if handles.add_noise.Value
                xmax = max(I2(:)).*2; %This *2 is to make up for the -.5 in the next var
                snrs = str2double(handles.noise_snr.String);
                noises = rand(size(I)).*xmax./snrs;
                I = I2+noises;
            end

            %%%% Calculate AE Signal
            for j = 1:length(ctime) %current time
                for i = 1:kgrid.Nt
                    if handles.uself.Value
                        for k = 1:size(L,3) %Lead Field
                            %                     V(:,:,i,j,k) = pressure(:,:,i).*(I(:,:,i,j).*L(:,:,k)).*med.k.*med.rho;
                            V2(:,:,i,j,k) = P(:,:,i).*I(:,:,j).*LFM(:,:,k).*med.k.*med.rho;
                            M(i,k,:) = squeeze(sum(sum(V(:,:,i,j,k))));
                        end
                    else
                        %V2(:,:,m,i,j) = P(:,:,m,i).*I(:,:,j).*med.k.*med.rho;
                        V2 = P(:,:,m,i).*I(:,:,j).*med.k.*med.rho;
                        M(i,m) = squeeze(sum(sum(V2)));
                    end
                end
            end
            multiWaitbar('Solving Forward Solution',m/size(P,3));
        end
        %  V2 = squeeze(V2);
        %         if handles.add_noise.Value
        %             xmax = max(M(:)).*2; %This *2 is to make up for the -.5 in the next var
        %             %    noises = (rand([size(V2,1),size(V2,2)])-.5).*xmax; % Too large
        %             snrs = str2double(handles.noise_snr.String);
        %             for i = 1:size(M,2)
        %
        %                 Mn(:,i) = M(:,i) +  (rand([size(M,1),1])-.5).*xmax./snrs;
        %
        %                 multiWaitbar('Adding Noise',i/size(V2,3));
        %             end
        %         else
        %             Mn = M;
        %         end
        %Xnoise = V2 + noise_var;
        %         V = squeeze(sum(V2,3));
         assignin('base','M',M);

         %%% Reconstruct Image
         if handles.reconstruct.Value
             V2 = zeros(size(P));
             for j = 1:size(P,3)
                 for i  = 1:size(P,4)

                     V2(:,:,j,i) = P(:,:,j,i) .* M(i,j);
                 end
                 multiWaitbar('Reconstructing Image',j/size(P,3));
             end

             % assignin('base','V2',V2);
             V = squeeze(sum(V2,3));
             V = squeeze(sum(V,3))';
             V2 = sum(V2,4);
             assignin('base','V4',V2);
         end
         multiWaitbar('CLOSEALL');
        %         if handles.add_noise.Value
        %             Mn = squeeze(mean(Mn,3));
        %             assignin('base','Mn',Mn);
        %         end
       
%        V = V2;
end

        % V2 = ones(size(V));
        % for i = 1:kgrid.Nt
        %     for j = 1:length(ctime)
        %         V2(:,:,i,j) = V(:,:,i,j)./pressure(:,:,i);
%
%     end
%     multiWaitbar('Solving Forward Solution',i/kgrid.Nt);
% end
% V = V2;
% for i = 1:kgrid.Nt
%     for j = 1:length(ctime)
%         for k = 1:size(I,1)
%         V(k,:,i,j) = pressure(k,:,i).*I(k,:,i,j).*med.k.*med.rho;
%         end
%     end
%     multiWaitbar('Solving Forward Solution',i/kgrid.Nt);
% end
if handles.reconstruct.Value
    assignin('base','V',V');
end
% assignin('base','M',M);
multiWaitbar('CLOSEALL');
% assignin('base','I',I);
% figure;
% a = max(V(:));
% b = min(V(:));
% for i = 1:kgrid.Nt
%     imagesc(V(:,:,i),[b a])
%     pause(0.05);
% end

% Localize V by considering considering t and c
% Export V as B-Mode

% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function medk_Callback(hObject, eventdata, handles)
% hObject    handle to medk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of medk as text
%        str2double(get(hObject,'String')) returns contents of medk as a double


% --- Executes during object creation, after setting all properties.
function medk_CreateFcn(hObject, eventdata, handles)
% hObject    handle to medk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function medrho_Callback(hObject, eventdata, handles)
% hObject    handle to medrho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of medrho as text
%        str2double(get(hObject,'String')) returns contents of medrho as a double


% --- Executes during object creation, after setting all properties.
function medrho_CreateFcn(hObject, eventdata, handles)
% hObject    handle to medrho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eleck_Callback(hObject, eventdata, handles)
% hObject    handle to eleck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eleck as text
%        str2double(get(hObject,'String')) returns contents of eleck as a double


% --- Executes during object creation, after setting all properties.
function eleck_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eleck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function elecrho_Callback(hObject, eventdata, handles)
% hObject    handle to elecrho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of elecrho as text
%        str2double(get(hObject,'String')) returns contents of elecrho as a double


% --- Executes during object creation, after setting all properties.
function elecrho_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elecrho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function elecx_Callback(hObject, eventdata, handles)
% hObject    handle to elecx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of elecx as text
%        str2double(get(hObject,'String')) returns contents of elecx as a double


% --- Executes during object creation, after setting all properties.
function elecx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elecx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function elecy_Callback(hObject, eventdata, handles)
% hObject    handle to elecy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of elecy as text
%        str2double(get(hObject,'String')) returns contents of elecy as a double


% --- Executes during object creation, after setting all properties.
function elecy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elecy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function medx_Callback(hObject, eventdata, handles)
% hObject    handle to medx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of medx as text
%        str2double(get(hObject,'String')) returns contents of medx as a double


% --- Executes during object creation, after setting all properties.
function medx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to medx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function medy_Callback(hObject, eventdata, handles)
% hObject    handle to medy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of medy as text
%        str2double(get(hObject,'String')) returns contents of medy as a double


% --- Executes during object creation, after setting all properties.
function medy_CreateFcn(hObject, eventdata, ~)
% hObject    handle to medy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
P = evalin('base','pressure');
kgrid = evalin('base','kgrid');
%P = permute(P,[2,1,3,4]);
xsize = size(P,1);
ysize = size(P,2);
t = size(P,3);
yax = linspace(0,kgrid.dx*xsize,xsize).*1000;
xax = linspace(0,kgrid.dy*ysize,ysize).*1000;
tax = linspace(0,kgrid.dt*t,t).*100000;

xrange = str2num(handles.xrange.String);
yrange = str2num(handles.yrange.String);
trange = str2num(handles.trange.String);
xpt = str2num(handles.xpoint.String);
ypt = str2num(handles.ypoint.String);
tpt = str2num(handles.tpoint.String);

if handles.ef.Value
    figure(9)
else
    axes(handles.axes1);
end
if ~handles.moviebox.Value
    
    P1 = P(xrange(1):xrange(2),yrange(1):yrange(2),tpt(1));
    if handles.show_envelope.Value
        %       P1 = P1.^2./(1.48*1e6);
        P1 = abs(hilbert(P1'))';
        %       P1 = envelope(P1);
        %       P1 = imgaussfilt(P1,1);
    end
    if isempty(handles.ae_clims.String)
        if handles.show_envelope.Value
            clims = [0 max(P1(:))];
        else
            clims = [-max(abs(P1(:))) max(abs(P1(:)))];
        end
    else
        clims = str2num(handles.ae_clims.String);
    end
    cmap = handles.colorm.String{handles.colorm.Value};
    Plot3(hObject, eventdata, handles, xax(yrange(1):yrange(2)), yax(xrange(1):xrange(2))-max(yax)/2, P1, clims, cmap)
% imagesc(xax(xrange(1):xrange(2)),yax(yrange(1):yrange(2)),P1,'ButtonDownFcn',{@PlotMetrics,handles})
% caxis([min(P1(:)) max(P1(:))]);
% xlabel('Depth mm')
% ylabel('Lateral mm')
% colormap(ones(64,3));
% colormap('gray');
% colorbar;
end
P2 = squeeze(P(xrange(1):xrange(2),ypt,trange(1):trange(2)));
P3 = squeeze(P(xpt,yrange(1):yrange(2),trange(1):trange(2)));

axes(handles.axes2)
imagesc(tax,xax,P2')
xlabel('X mm')
ylabel('T us')

axes(handles.axes3)
imagesc(tax,yax,P3')
xlabel('Y mm')
ylabel('T us')

if handles.moviebox.Value
     if handles.save_movie.Value
        clims = [-1*max(P(:))*.9 max(P(:))*.9];
        file = handles.record_filename.String;
        path = handles.record_dir.String;
        v = VideoWriter(fullfile(path,file));
        v.FrameRate = str2double(handles.sim_fps.String);
        open(v);
        p = 'n';
        figure(100)
        T = tpt(1):tpt(2);
        for i = 1:length(T)
            imagesc(P(:,:,T(i)),'ButtonDownFcn',{@PlotMetrics,handles});
            title([p ' = ' num2str(T(i)*kgrid.dt*1e6) ' \mus']);
            colormap('gray');
            text(1,15,[p ' = ' num2str(T(i)*kgrid.dt*1e6) ' \mus'],'Color','white','FontSize',14);  
            caxis(clims)
            xlabel('Depth mm')
            ylabel('Lateral mm')
            frame = getframe;
            writeVideo(v,frame);
        end
        close(v)
    else
    axes(handles.axes1);
    xlabel('Depth mm')
    ylabel('Lateral mm')
    T = tpt(1):tpt(2);
    imagesc(xax,yax,P(xrange(1):xrange(2),yrange(1):yrange(2),T(1)),[min(P(:)) max(P(:))]);
    for i = 2:length(T)
%       imagesc(xax,yax,P(xrange(1):xrange(2),yrange(1):yrange(2),T(i)));
        handles.axes1.Children.CData = P(xrange(1):xrange(2),yrange(1):yrange(2),T(i));
        drawnow;
        title(num2str(i))
        pause(str2double(handles.movie_pause.String));
    end
     xlabel('Depth mm')
    ylabel('Lateral mm')
     end
end
        

    
    
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
P = evalin('base','I');
kgrid = evalin('base','kgrid');
xsize = size(P,1);
ysize = size(P,2);
t = size(P,3);
yax = linspace(0,kgrid.dx*xsize,xsize).*1000;
xax = linspace(0,kgrid.dy*ysize,ysize).*1000;
c = size(P,4);
cax = linspace(0,str2num(handles.cdt.String),c);

xrange = str2num(handles.xrange.String);
yrange = str2num(handles.yrange.String);
trange = str2num(handles.trange.String);
crange = str2num(handles.crange.String);

xpt = str2num(handles.xpoint.String);
ypt = str2num(handles.ypoint.String);
tpt = str2num(handles.tpoint.String);
cpt = str2num(handles.cpoint.String);

if handles.ef.Value
    figure(11);
else
axes(handles.axes1);
end
if ~handles.moviebox.Value
    P1 = P(xrange(1):xrange(2),yrange(1):yrange(2),1,cpt(1));
    imagesc(xax(yrange(1):yrange(2)),yax(xrange(1):xrange(2)),P1,'ButtonDownFcn',{@PlotMetrics,handles});
    caxis([min(P(:)) max(P(:))])
    colormap(ones(64,3));
colormap('hot');
    colorbar;
     xlabel('Depth mm')
    ylabel('Lateral mm')
end
if handles.slowtimebox.Value
    tax = cax;
    P2 = squeeze(P(xrange(1):xrange(2),ypt(1),tpt(1),crange(1):crange(2)));
    P3 = squeeze(P(xpt(1),yrange(1):yrange(2),tpt(1),crange(1):crange(2)));
    colorbar;
else
%     P2 = squeeze(P(xrange(1):xrange(2),ypt,trange(1):trange(2),cpt(1)));
%     P3 = squeeze(P(xpt,yrange(1):yrange(2),trange(1):trange(2),cpt(1)));
%     colorbar;
end


% axes(handles.axes2)
% imagesc(tax,xax,P2)
%   xlabel('T us')
%     ylabel('X mm')
% 
% axes(handles.axes3)
% imagesc(tax,yax,P3)
%   xlabel('T us')
%     ylabel('Y mm')

if handles.moviebox.Value
        if handles.save_movie.Value
        clims = [min(P(:)) max(P(:))];
        file = handles.record_filename.String;
        path = handles.record_dir.String;
        v = VideoWriter(fullfile(path,file));
        v.FrameRate = str2double(handles.sim_fps.String);
        open(v);
        p = 'n';
        figure(100)
        T = tpt(1):tpt(2);
        for i = 1:length(T)
            imagesc(P(:,:,T(i)));
            title([p ' = ' num2str(T(i)*kgrid.dt*1e6) ' \mus']);
            colormap('hot');
            text(1,15,[p ' = ' num2str(T(i)*kgrid.dt*1e6) ' \mus'],'Color','white','FontSize',14);  
            caxis(clims)
            xlabel('Depth mm')
            ylabel('Lateral mm')
            frame = getframe;
            writeVideo(v,frame);
        end
        close(v)
    else
    clims = [min(P(:)) max(P(:))];
    axes(handles.axes1);
    if handles.slowtimebox.Value
        T = cpt(1):cpt(2);
        P1 = squeeze(P(xrange(1):xrange(2),yrange(1):yrange(2),tpt(1),crange(1):crange(2)));
    else
    T = tpt(1):tpt(2);
    P1 = squeeze(P(xrange(1):xrange(2),yrange(1):yrange(2),trange(1):trange(2),cpt(1)));
    end
    for i = 1:length(T)
        title(num2str(i));
        imagesc(xax,yax,P1(:,:,T(i)),clims);
        colorbar;
        pause(str2double(handles.movie_pause.String));
    end
    xlabel('Depth mm')
    ylabel('Lateral mm')
        end
end

% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
ptype = handles.pressuremenu.String{handles.pressuremenu.Value};
switch ptype
    case  {'sensor_data','pressure','Hadamard','planes'}
        P = evalin('base','V');
        L = handles.LF_Menu.Value;
        if ndims(P) == 4
            P = P(:,:,:,L);
        elseif ndims(P) == 5
            P = P(:,:,:,:,L);
        end
        kgrid = evalin('base','kgrid');
        xsize = size(P,1);
        ysize = size(P,2);
        t = size(P,3);
        yax = linspace(0,kgrid.dx*xsize,xsize).*1000;
        xax = linspace(0,kgrid.dy*ysize,ysize).*1000;
        c = size(P,4);
        tax = linspace(0,kgrid.dt*t,t)*1000000;
        cax = linspace(0,str2num(handles.cdt.String),c);
        
        xrange = str2num(handles.xrange.String);
        yrange = str2num(handles.yrange.String);
        trange = str2num(handles.trange.String);
        crange = str2num(handles.crange.String);
        
        xpt = str2num(handles.xpoint.String);
        ypt = str2num(handles.ypoint.String);
        tpt = str2num(handles.tpoint.String);
        cpt = str2num(handles.cpoint.String);
        
        if handles.show_envelope.Value
            if handles.moviebox.Value
                for i = 1:size(P,3)
                    P(:,:,i) = abs(hilbert(P(:,:,i)'))';
                    if isempty(handles.ae_clims.String)
                        clims = [min(P(:)) max(P(:))];
                    else
                        clims = str2num(handles.ae_clims.String);
                    end
                    
                end
            else
                P(:,:,tpt) = abs(hilbert(P(:,:,tpt)'))';
                if isempty(handles.ae_clims.String)
                    clims = [min(min(P(:,:,tpt))) max(max(P(:,:,tpt)))];
                else
                    clims = str2num(handles.ae_clims.String);
                end
            end
            if strcmp(handles.colorm.String{handles.colorm.Value},'default')
                cmaps = 'hot';
            else
                cmaps = handles.colorm.String{handles.colorm.Value};
            end
            %     clims = [min(P(:)) max(P(:))];
        else
            if strcmp(handles.colorm.String{handles.colorm.Value},'default')
                cmaps = 'hotcold';
            else
                cmaps = handles.colorm.String{handles.colorm.Value};
            end
            if isempty(handles.ae_clims.String)
                clims = [-1*max(abs(P(:))) max(abs(P(:)))];
            else
                clims = str2num(handles.ae_clims.String);
            end
            
        end
        
        
        if handles.ef.Value
            figure(12);
        else
            axes(handles.axes1);
        end
        if ~handles.moviebox.Value
            P1 = P(xrange(1):xrange(2),yrange(1):yrange(2),tpt(1),cpt(1));
            %     imagesc(xax(xrange(1):xrange(2)),yax(yrange(1):yrange(2)),P1,'ButtonDownFcn',{@PlotMetrics,handles});
            Plot3(hObject, eventdata, handles, xax(yrange(1):yrange(2)), yax(xrange(1):xrange(2)), P1, clims, cmaps)
            %     caxis(clims);
            %     colormap(ones(64,3));
            %     colormap(cmaps);
            %     colorbar;
            %     xlabel('Depth mm')
            %     ylabel('Lateral mm')
        end
        if handles.slowtimebox.Value
            tax = cax;
            P2 = squeeze(P(xrange(1):xrange(2),ypt(1),tpt(1),crange(1):crange(2)));
            P3 = squeeze(P(xpt(1),yrange(1):yrange(2),tpt(1),crange(1):crange(2)));
            colorbar;
        else
            P2 = squeeze(P(xrange(1):xrange(2),ypt(1),trange(1):trange(2),cpt(1)));
            P3 = squeeze(P(xpt(1),yrange(1):yrange(2),trange(1):trange(2),cpt(1)));
            colorbar;
        end
        
        
%         axes(handles.axes2)
%         imagesc(tax,xax,P2)
%         xlabel('T us')
%         ylabel('X mm')
%         
%         axes(handles.axes3)
%         imagesc(tax,yax,P3)
%         xlabel('T us')
%         ylabel('Y mm')
        
        if handles.moviebox.Value
            if handles.save_movie.Value
                clims = [-1*max(P(:)) max(P(:))];
                file = handles.record_filename.String;
                path = handles.record_dir.String;
                v = VideoWriter(fullfile(path,file));
                v.FrameRate = str2double(handles.sim_fps.String);
                open(v);
                p = 'n';
                figure(100)
                T = tpt(1):tpt(2);
                for i = 1:length(T)
                    imagesc(P(:,:,T(i)));
                    title([p ' = ' num2str(T(i)*kgrid.dt*1e6) ' \mus']);
                    colormap('hotcold');
                    text(1,15,[p ' = ' num2str(T(i)*kgrid.dt*1e6) ' \mus'],'Color','white','FontSize',14);
                    caxis(clims)
                    xlabel('Depth mm')
                    ylabel('Lateral mm')
                    frame = getframe;
                    writeVideo(v,frame);
                end
                close(v)
            else
                clims = [-1*max(P(:)) max(P(:))];
                axes(handles.axes1);
                if handles.slowtimebox.Value
                    T = cpt(1):cpt(2);
                    P1 = squeeze(P(xrange(1):xrange(2),yrange(1):yrange(2),tpt,crange(1):crange(2)));
                else
                    T = tpt(1):tpt(2);
                    P1 = squeeze(P(xrange(1):xrange(2),yrange(1):yrange(2),trange(1):trange(2),cpt));
                end
                %         P4 = ones(size(P1),'gpuArray');
                %         P1 = P4.*P1;
                tic
                for i = 1:length(T)
                    imagesc(xax,yax,P1(:,:,T(i)),clims);
                    colorbar;
                    title(num2str(i));
                    pause(str2double(handles.movie_pause.String));
                end
                toc
                xlabel('Depth mm')
                ylabel('Lateral mm')
            end
        end
    case {'convolution','multi_conv'}
        P = evalin('base','M2');
        L = handles.LF_Menu.Value;
        
        P = P(:,:,L);
        
        kgrid = evalin('base','kgrid');
        xsize = size(P,2);
        ysize = size(P,1);
        if handles.show_envelope.Value
            
            P = abs(hilbert(P));
            if isempty(handles.ae_clims.String)
            clims = [min(min(P)) max(max(P))];
            else 
                clims = str2num(handles.ae_clims.String);
            end
            
            if strcmp(handles.colorm.String{handles.colorm.Value},'default')
                cmaps = 'hot';
            else
                cmaps = handles.colorm.String{handles.colorm.Value};
            end
            %     clims = [min(P(:)) max(P(:))];
        else
            if strcmp(handles.colorm.String{handles.colorm.Value},'default')
                cmaps = 'hotcold';
            else
                cmaps = handles.colorm.String{handles.colorm.Value};
            end
            if isempty(handles.ae_clims.String)
                clims = [-1*max(abs(P(:))) max(abs(P(:)))];
            else
                clims = str2num(handles.ae_clims.String);
            end
            
        end
        
        
        if handles.ef.Value
            figure(12);
        else
            axes(handles.axes1);
        end
        if ~handles.moviebox.Value
            ysize = (size(P,2)-1)/2;
            kgrid = evalin('base','kgrid');
            xsize = kgrid.Nt*kgrid.dt*1e6;
%             imagesc(linspace(0,xsize,size(P,2)),linspace(ysize*-1,ysize,size(P,2)),P');
            Plot3(hObject, eventdata, handles, linspace(0,kgrid.dt*1e6*1.480*size(P,1),size(P,1)),linspace(-4.5,4.5,30),P',clims,cmaps);
            %     P1 = P(xrange(1):xrange(2),yrange(1):yrange(2),tpt(1),cpt(1));
            %     %     imagesc(xax(xrange(1):xrange(2)),yax(yrange(1):yrange(2)),P1,'ButtonDownFcn',{@PlotMetrics,handles});
            %     Plot3(hObject, eventdata, handles, xax(yrange(1):yrange(2)), yax(xrange(1):xrange(2)), P1, clims, cmaps)
            %     caxis(clims);
            %     colormap(ones(64,3));
            %     colormap(cmaps);
            %     colorbar;
            %     xlabel('Depth mm')
            %     ylabel('Lateral mm')
        end
end



        

function Plot3(hObject, eventdata, handles, xax, yax, P, clims, cmaps)
if handles.db_box.Value
    if isempty(handles.dbmax.String)
    P = real(20*log10(P./max(P(:))));
    else
         P = real(20*log10(P./str2double(handles.dbmax.String)));
    end
end
imagesc(xax,yax,P,'ButtonDownFcn',{@PlotMetrics,handles});
if clims(1) < clims(2) 
caxis(clims);
end
 colormap(ones(64,3));
 colormap(cmaps);
 colorbar;
 xlabel('Depth mm')
 ylabel('Lateral mm')
 
     



function PlotMetrics(hObject, eventdata, handles)
pt = get(gca,'currentpoint');
x = pt(1,1);
y = pt(1,2);
if handles.exact.Value
    if ~isempty(handles.exactx.String)
        x = str2double(handles.exactx.String);
    end
    if ~isempty(handles.exacty.String)
        y = str2double(handles.exacty.String);
    end
end
xrange = str2double(handles.metrics_yr.String)/2;
yrange = str2double(handles.metrics_xr.String)/2;
Xlims(1) = find(hObject.XData >= x-xrange,1);
Xlims(2) = find(hObject.XData <=x+xrange,1,'last');
Ylims(1) = find(hObject.YData >= y-yrange,1);
Ylims(2) = find(hObject.YData <=y+yrange,1,'last');
psource = handles.pressuremenu.String{handles.pressuremenu.Value};


    p = hObject.CData;
P = p(Ylims(1):Ylims(2),Xlims(1):Xlims(2));
pmax = find(P == max(P(:)),1);
if handles.exact.Value
    kgrid = evalin('base','kgrid');
    y_vec = (kgrid.x_vec)*1e3;
     x_vec = (kgrid.y_vec - min(kgrid.y_vec))*1e3;
     x_vec = x_vec(Xlims(1):Xlims(2));
     y_vec = y_vec(Ylims(1):Ylims(2));
    xval = find(abs(y-y_vec) == min(abs(y-y_vec)),1);
    yval = find(abs(x-x_vec) == min(abs(x-x_vec)),1);
%     xval = 100;
%     yval = 51;
else
    yval = ceil(pmax/size(P,1));
    xval = mod(pmax,size(P,1));
end
 figure(19);
plot(hObject.XData(1:size(P,2)),P(xval,:),'k','LineWidth',2);
title('Depth')
figure(20);
plot(hObject.YData(1:size(P,1)),P(:,yval),'k','LineWidth',2);
title('Lateral')

F = double(P(:,yval));
F = resample(F,10,1);
fmax = max(F);
fhalf = fmax/2;
fmid = find(F == fmax,1);
f1 = find(F(1:fmid) <= fhalf,1,'last');
f2 = find(F(fmid:end) <= fhalf,1);
f2 = f2+fmid;

G = double(P(xval,:));
G = resample(G,10,1);
gmax = max(G);
ghalf = gmax/2;
gmid = find(G == gmax,1);
g1 = find(G(1:gmid) <= ghalf,1,'last');
g2 = find(G(gmid:end) <= ghalf,1);
g2 = g2+gmid;


switch psource
    case 'convolution'
        spacing = evalin('base','spacing');
fwhm = (f2-f1)*spacing;
gwhm = g2-g1;
    case {'sensor_data','pressure'}
        fwhm = f2-f1;
        gwhm = g2-g1;
end
kgrid = evalin('base','kgrid');
gwhm = gwhm*kgrid.dy*100;
fwhm = fwhm*kgrid.dx*100;
set(handles.lateraltext,'String',fwhm);
set(handles.axialtext,'String',gwhm);
set(handles.peaktext,'String',max(P(:)));
disp(['Lateral FWHM is ' num2str(fwhm) ' mm']);
disp(['Axial FWHM is ' num2str(gwhm) ' mm']);
disp(['Peak value = ' num2str(max(P(:)))]);
int = sum(P(:));
set(handles.integraltext,'String',int);
disp(['Integral = ' num2str(int)]);
if handles.comp_box.Value
    comp_T = evalin('base','comp_T');
    Pnorm = P./(max(P(:)));
    comp_Tnorm = comp_T./max(comp_T(:));
    
    ME1 = abs(Pnorm(1:size(comp_Tnorm),:)-comp_Tnorm);
    ME = sum(ME1(:));
    set(handles.meanerrortext,'String',ME);
    
    disp(['Mean Error = ' num2str(ME)]);
end

if handles.ef.Value
    figure(21)
else
axes(handles.axes2)
end
yloc = abs(hObject.YData-y);
yloc = find(yloc == min(yloc));
plot(hObject.XData,hObject.CData(yloc,:),'k','LineWidth',2);

if handles.ef.Value
    figure(22)
else
axes(handles.axes3)
hold off
cla reset
end
xloc = abs(hObject.XData-x);
xloc = find(xloc == min(xloc));
switch psource
    case 'convolution'
        hObject.YData = hObject.YData./(kgrid.Nx/size(p,1)).*spacing;
end
if handles.template_box.Value
    assignin('base','comp_T',P);
    set(handles.template_box,'Value',0);
end
plot(hObject.YData,hObject.CData(:,xloc),'k','LineWidth',2);
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function xrange_Callback(hObject, eventdata, handles)
% hObject    handle to xrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xrange as text
%        str2double(get(hObject,'String')) returns contents of xrange as a double


% --- Executes during object creation, after setting all properties.
function xrange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yrange_Callback(hObject, eventdata, handles)
% hObject    handle to yrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yrange as text
%        str2double(get(hObject,'String')) returns contents of yrange as a double


% --- Executes during object creation, after setting all properties.
function yrange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trange_Callback(hObject, eventdata, handles)
% hObject    handle to trange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trange as text
%        str2double(get(hObject,'String')) returns contents of trange as a double


% --- Executes during object creation, after setting all properties.
function trange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xpoint_Callback(hObject, eventdata, handles)
% hObject    handle to xpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xpoint as text
%        str2double(get(hObject,'String')) returns contents of xpoint as a double


% --- Executes during object creation, after setting all properties.
function xpoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ypoint_Callback(hObject, eventdata, handles)
% hObject    handle to ypoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ypoint as text
%        str2double(get(hObject,'String')) returns contents of ypoint as a double


% --- Executes during object creation, after setting all properties.
function ypoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ypoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tpoint_Callback(hObject, eventdata, handles)
% hObject    handle to tpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tpoint as text
%        str2double(get(hObject,'String')) returns contents of tpoint as a double


% --- Executes during object creation, after setting all properties.
function tpoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in holdmax.
function holdmax_Callback(hObject, eventdata, handles)
% hObject    handle to holdmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of holdmax


% --- Executes on button press in max_global.
function max_global_Callback(hObject, eventdata, handles)
% hObject    handle to max_global (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of max_global


% --- Executes on button press in twodsensor.
function twodsensor_Callback(hObject, eventdata, handles)
% hObject    handle to twodsensor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of twodsensor



function camp_Callback(hObject, eventdata, handles)
% hObject    handle to camp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of camp as text
%        str2double(get(hObject,'String')) returns contents of camp as a double


% --- Executes during object creation, after setting all properties.
function camp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to camp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cfreq_Callback(hObject, eventdata, handles)
% hObject    handle to cfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cfreq as text
%        str2double(get(hObject,'String')) returns contents of cfreq as a double


% --- Executes during object creation, after setting all properties.
function cfreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cdelay_Callback(hObject, eventdata, handles)
% hObject    handle to cdelay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cdelay as text
%        str2double(get(hObject,'String')) returns contents of cdelay as a double


% --- Executes during object creation, after setting all properties.
function cdelay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cdelay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cdur_Callback(hObject, eventdata, handles)
% hObject    handle to cdur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cdur as text
%        str2double(get(hObject,'String')) returns contents of cdur as a double


% --- Executes during object creation, after setting all properties.
function cdur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cdur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in TRbutton.
function TRbutton_Callback(hObject, eventdata, handles)
if handles.tr_all.Value
    %     TR = evalin('base','rev');
    %     col = 'krgbym';
    %     T2 = TR-min(TR);
    %     if handles.ef.Value
    %         for i = 1:size(TR,2)
    %             figure(30+i)
    % %             T2 = TR(:,i)-min(TR(:,i));
    %             plot(T2(:,i),col(i),'LineWidth',2.5)
    %             xlabel('Element Number')
    %             ylabel('Delay Added [\mus]')
    %             ylim([0 ceil(max(T2(:)))])
    %             title(num2str(i))
    %         end
    %     end
     if handles.tr_apod.Value
         TR = evalin('base','Apods');
          source = evalin('base','source');
        kgrid = evalin('base','kgrid');
        aM = TR/min(TR(:));
        iMax = sum(abs(source.p(:)));
        for i = 1:size(source.p,1)
            source.p(i,:) = source.p(i,:)./aM(i);
        end
        
        %Normalize to 1
        newMax = sum(abs(source.p(:)));
        aMax = TR/max(TR(:));
        aMean = 1/mean(aMax);
        %         source.p = source.p.*aMean;
        source.p = source.p.*(iMax/newMax);
        assignin('base','source',source);
     else
         TR = evalin('base','TR2');
         source = evalin('base','source');
         kgrid = evalin('base','kgrid');
         padN = round((max(TR)-min(TR))/(kgrid.dt*1e6));
         
         for i = 1:size(source.p,1)
             shiftN = round((TR(i)-min(TR))/(kgrid.dt*1e6));
             p2 = padarray(source.p(i,:)',padN,'pre');
             p3(i,:) = circshift(p2,-shiftN);
         end
         source.p = p3;
         assignin('base','source',source);
         if handles.ef.Value
             figure(8)
             T2 = TR-min(TR);
             plot(T2,'k','LineWidth',2.5)
             xlabel('Element Number')
             ylabel('Delay Added [\mus]')
         end
     end
         
else
    if handles.tr_apod.Value
        R = inputdlg('Input apods');
        R2 = R{1};
        if isempty(R2)
            TR = evalin('base','Apods');
        else
            TR = str2num(R2);
        end
        source = evalin('base','source');
        kgrid = evalin('base','kgrid');
        aM = TR/min(TR(:));
        iMax = sum(abs(source.p(:)));
        for i = 1:size(source.p,1)
            source.p(i,:) = source.p(i,:)./aM(i);
        end
        
        %Normalize to 1
        newMax = sum(abs(source.p(:)));
        aMax = TR/max(TR(:));
        aMean = 1/mean(aMax);
        %         source.p = source.p.*aMean;
        source.p = source.p.*(iMax/newMax);
        assignin('base','source',source);
    else
        R = inputdlg('Input delays');
        R2 = R{1};
        if isempty(R2)
            TR = evalin('base','TR');
        else
            TR = str2num(R2);
        end
        source = evalin('base','source');
        kgrid = evalin('base','kgrid');
        padN = round((max(TR)-min(TR))/(kgrid.dt*1e6));
        
        for i = 1:size(source.p,1)
            shiftN = round((TR(i)-min(TR))/(kgrid.dt*1e6));
            p2 = padarray(source.p(i,:)',padN,'pre');
            p3(i,:) = circshift(p2,-shiftN);
        end
        source.p = p3;
        assignin('base','source',source);
        if handles.ef.Value
            figure(8)
            T2 = TR-min(TR);
            plot(T2,'k','LineWidth',2.5)
            xlabel('Element Number')
            ylabel('Delay Added [\mus]')
        end
    end
end


% hObject    handle to TRbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in moviebox.
function moviebox_Callback(hObject, eventdata, handles)
% hObject    handle to moviebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of moviebox


% --- Executes on button press in ae_max.
function ae_max_Callback(hObject, eventdata, handles)
P = evalin('base','pressure');
x = size(P,1);
y = size(P,2);
t = size(P,3);
c = size(P,4);
set(handles.xrange,'String',num2str([1 x]));
set(handles.yrange,'String',num2str([1 y]));
set(handles.trange,'String',num2str([1 t]));
set(handles.crange,'String',num2str([1 c]));


% hObject    handle to ae_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function movie_pause_Callback(hObject, eventdata, handles)
% hObject    handle to movie_pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of movie_pause as text
%        str2double(get(hObject,'String')) returns contents of movie_pause as a double


% --- Executes during object creation, after setting all properties.
function movie_pause_CreateFcn(hObject, eventdata, handles)
% hObject    handle to movie_pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cdt_Callback(hObject, eventdata, handles)
% hObject    handle to cdt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cdt as text
%        str2double(get(hObject,'String')) returns contents of cdt as a double


% --- Executes during object creation, after setting all properties.
function cdt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cdt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function crange_Callback(hObject, eventdata, handles)
% hObject    handle to crange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of crange as text
%        str2double(get(hObject,'String')) returns contents of crange as a double


% --- Executes during object creation, after setting all properties.
function crange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to crange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cpoint_Callback(hObject, eventdata, handles)
% hObject    handle to cpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cpoint as text
%        str2double(get(hObject,'String')) returns contents of cpoint as a double


% --- Executes during object creation, after setting all properties.
function cpoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in slowtimebox.
function slowtimebox_Callback(hObject, eventdata, handles)
% hObject    handle to slowtimebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of slowtimebox


% --- Executes on button press in checkbox19.
function checkbox19_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox19


% --- Executes on selection change in electroidN.
function electroidN_Callback(hObject, eventdata, handles)
% hObject    handle to electroidN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns electroidN contents as cell array
%        contents{get(hObject,'Value')} returns selected item from electroidN


% --- Executes during object creation, after setting all properties.
function electroidN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to electroidN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function make_electroid_Callback(hObject, eventdata, handles)
% hObject    handle to make_electroid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of make_electroid as text
%        str2double(get(hObject,'String')) returns contents of make_electroid as a double


% --- Executes during object creation, after setting all properties.
function make_electroid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to make_electroid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Elec_Button.
function Elec_Button_Callback(hObject, eventdata, handles)
Electroids;
% hObject    handle to Elec_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function Elec_Button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Elec_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in use_elec.
function use_elec_Callback(hObject, eventdata, handles)
% hObject    handle to use_elec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_elec


% --- Executes on button press in autoelecmedia.
function autoelecmedia_Callback(hObject, eventdata, handles)
kgrid = evalin('base','kgrid');
y = kgrid.Ny;
x = kgrid.Nx;
set(handles.medx,'String',x);
set(handles.medy,'String',y);
p = evalin('base','sensor_data.p');
pressure = reshape(p,[x y size(p,2)]); 
assignin('base','pressure',pressure);


% hObject    handle to autoelecmedia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in show_envelope.
function show_envelope_Callback(hObject, eventdata, handles)
% hObject    handle to show_envelope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of show_envelope


% --- Executes on button press in gauss_apod.
function gauss_apod_Callback(hObject, eventdata, handles)
% hObject    handle to gauss_apod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gauss_apod


% --- Executes on button press in plotFWHM.
function plotFWHM_Callback(hObject, eventdata, handles)
% [x,y] = ginput(1);

% x = round(x);
% data = sensor_data.p;
% ymax = data(x);
% ymin =
% yhalf = ymax/2;
if handles.twodsensor.Value
    M = evalin('base','M2');
    M = M(:,:,handles.LF_Menu.Value);
    if handles.db_box.Value
        if isempty(handles.dbmax.String)
            P = real(20*log10(M./max(M(:))));
        else
            P = real(20*log10(P./str2double(handles.dbmax.String)));
        end
    end
    kgrid = evalin('base','kgrid');
    sig = find(P >= -6);
    xsize = 0.3;
    ysize = kgrid.dt*1e6*1.485;
    area = length(sig)*xsize*ysize;
    peak = max(abs(M(:)));
    disp(['Area of signal is ' num2str(area) ' with Peak ' num2str(peak)]);
else

    sensor_data = evalin('base','sensor_data');
    kgrid = evalin('base','kgrid');
    snum = handles.sensormenu.Value;
    t = str2double(handles.tp.String);
    p = sensor_data.p(:,t);

    p = reshape(p,[sqrt(size(p,1)),sqrt(size(p,1)),size(p,2)]);

    p = envelope(p);

    ymax = max(p(:));
    yhalf = ymax/2;
    F = find(p == ymax,1);
    col = ceil(F/size(p,1));
    row = mod(F,size(p,1));
    axes(handles.axes1)
    imagesc(p);
    pvert = p(:,col);
    phorz = p(row,:);
    pvert = resample(pvert,10,1);
    axes(handles.axes2)
    plot(pvert);
    phorz = resample(phorz,10,1);
    axes(handles.axes3)
    hold off
    cla reset
    plot(phorz);
    s2 = snum*10;
    v_first = pvert(1:row*10);
    v_second = pvert(row*10+1:end);
    h_first = phorz(1:col*10);
    h_second = phorz(col*10+1:end);
    x1 = find(v_first <= yhalf,1,'last');
    x2 = find(v_second <= yhalf,1)+row*10;
    kgrid = evalin('base','kgrid');
    dy = kgrid.dy/10;
    fwhm = (x2-x1)*dy*1e3;
    dx = kgrid.dx/10;
    y1 = find(h_first <= yhalf,1,'last');
    y2 = find(h_second <= yhalf,1)+col*10;
    fwhm2 = (y2-y1)*dx*1e3;
    set(handles.peaktext,'String',ymax);
    set(handles.lateraltext,'String',fwhm2);
    set(handles.axialtext,'String',fwhm);
    disp(['FWHM is ' num2str(fwhm) ' mm in depth, ' num2str(fwhm2) ' mm laterally, and peak is ' num2str(ymax)]);
end

% hObject    handle to plotFWHM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in find_max_t.
function find_max_t_Callback(hObject, eventdata, handles)
%sensor_data = evalin('base','sensor_data');
snum = handles.sensormenu.Value;
pressure = evalin('base','pressure');
kgrid = evalin('base','kgrid');
if handles.show_envelope.Value
    p = zeros(size(pressure));
        for i = 1:size(pressure,3)
        p(:,:,i) = abs(hilbert(pressure(:,:,i)'))';
        multiWaitbar('Enveloping',i/size(pressure,3));
        end
        multiWaitbar('CLOSEALL');
%     p = abs(pressure);
    [x y] = ind2sub([kgrid.Nx kgrid.Ny],snum);
    p2 = squeeze(p(x,y,:));
    t = find(p2 == max(p2),1);
else
    if ndims(pressure) == 3
        if handles.box3d.Value
            plane = handles.p_menu.String{handles.p_menu.Value};
            switch plane
                case 'XY'
                    [x y] = ind2sub([kgrid.Ny kgrid.Nz],snum);  
                case 'YZ'
                    [x y] = ind2sub([kgrid.Ny kgrid.Nz],snum);         
                case 'XZ'
                    [x y] = ind2sub([kgrid.Nz kgrid.Nx],snum);
            end
        else
            [x y] = ind2sub([kgrid.Nx kgrid.Ny],snum);
          
        end
         p = abs(squeeze(pressure(x,y,:)));
    else
        [x y z] = ind2sub([kgrid.Nx kgrid.Ny kgrid.Nz],snum);
        p = abs(squeeze(pressure(x,y,z,:)));
    end
    t = find(p == max(p),1);
end


kgrid = evalin('base','kgrid');
t2 = t*kgrid.dt*1e6;
disp(['t = ' num2str(t2) ' us']);

set(handles.tp,'String',t);
set(handles.tslide,'Value',t);
set(handles.tpoint,'String',t);
% show_tp_Callback(hObject, eventdata, handles)
% hObject    handle to find_max_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function tslide_Callback(hObject, eventdata, handles)
t = round(hObject.Value);
set(handles.tp,'String',t);
show_tp_Callback(hObject, eventdata, handles)

% hObject    handle to tslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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


% --- Executes on button press in autofill.
function autofill_Callback(hObject, eventdata, handles)
H = handles.automenu.String{handles.automenu.Value};
if ismember('kgrid',evalin('base','who'))
    kgrid = evalin('base','kgrid');
end
switch H
    case 'H235'
        set(handles.Nx,'String',300);
        set(handles.Ny,'String',240);
        set(handles.Nt,'String',600);
        set(handles.Nz,'String',100);
        set(handles.dx,'String',0.25);
        set(handles.dy,'String',0.27);
        set(handles.dz,'String',0.45)
        set(handles.dt_sim,'String',0.1);
        make_grid_Callback(hObject, eventdata, handles);
        set(handles.trans_freq,'String',0.6)
        set(handles.trans_cycles,'String',2)
        set(handles.gauss_apod,'Value',1)
        set(handles.wave_type,'Value',1)
        make_medium_Callback(hObject, eventdata, handles);
        set(handles.trans_x0,'String',10)
        set(handles.trans_y0,'String',120)
        set(handles.trans_z0,'String',50)
        set(handles.trans_focus,'String',num2str([150 120 50]))
        set(handles.box3d,'Value',1)
        set(handles.all_sensor,'Value',0)
        set(handles.s_plane,'String',50)
        make_sensor_Callback(hObject, eventdata, handles);
        set(handles.trans_elements,'String',126)
        set(handles.s_plane,'String',50)
        set(handles.az_ele,'String',18)
        set(handles.ela_ele,'String',7)
        set(handles.trans_active_elements,'String','1:126')
        set(handles.trans_x_kerf,'String',10)
        set(handles.trans_y_kerf,'String',10)
        set(handles.trans_length,'String',9)
        set(handles.trans_width,'String',9)
        set(handles.source_3d,'Value',1)
        set(handles.array_2d,'Value',1)
        set(handles.curve_el,'String',num2str([16 11 5 0 5 11 16]))
        set(handles.sim_3D,'Value',1)   
        set(handles.wave_type,'Value',2)
    case 'Hadamard2D'
        set(handles.Nx,'String',300);
        set(handles.Ny,'String',300);
        set(handles.Nt,'String',1800);
        set(handles.Nz,'String',0);
        set(handles.dx,'String',0.2);
        set(handles.dy,'String',0.2);
        set(handles.dz,'String',0.45)
        set(handles.dt_sim,'String',0.025);
        make_grid_Callback(hObject, eventdata, handles);
        set(handles.trans_freq,'String',2.5)
        set(handles.trans_cycles,'String',2)
        set(handles.gauss_apod,'Value',1)
        set(handles.wave_type,'Value',1)
     %   make_medium_Callback(hObject, eventdata, handles);
        set(handles.trans_x0,'String',150)
        set(handles.trans_y0,'String',10)
        set(handles.trans_z0,'String',50)
        set(handles.trans_focus,'String',num2str([150 120 50]))
        set(handles.box3d,'Value',0)
        set(handles.all_sensor,'Value',1)
        set(handles.s_plane,'String',0)
        make_sensor_Callback(hObject, eventdata, handles);
        set(handles.trans_elements,'String',16)
        set(handles.az_ele,'String',16)
        set(handles.ela_ele,'String',1)
        set(handles.trans_active_elements,'String','1:16')
        set(handles.trans_x_kerf,'String',9)
        set(handles.trans_y_kerf,'String',1)
        set(handles.trans_length,'String',1)
        set(handles.trans_width,'String',1)
        set(handles.source_3d,'Value',0)
        set(handles.array_2d,'Value',0)    
        set(handles.sim_3D,'Value',0)   
        set(handles.wave_type,'Value',1)
    case 'ReconKgrid6x6'
        set(handles.Nx,'String',400);
        set(handles.Ny,'String',400);
        set(handles.Nt,'String',1000);
        set(handles.Nz,'String',0);
        set(handles.dx,'String',0.15);
        set(handles.dy,'String',0.15);
        set(handles.dz,'String',0.15)
        set(handles.dt_sim,'String',0.05);
        make_grid_Callback(hObject, eventdata, handles);
    case 'ReconKgrid6x9'
               set(handles.Nx,'String',400);
        set(handles.Ny,'String',600);
        set(handles.Nt,'String',1400);
        set(handles.Nz,'String',0);
        set(handles.dx,'String',0.15);
        set(handles.dy,'String',0.15);
        set(handles.dz,'String',0.15)
        set(handles.dt_sim,'String',0.05);
        make_grid_Callback(hObject, eventdata, handles);
    case 'ReconKgrid5x4x8'
        set(handles.Nx,'String',320);
        set(handles.Ny,'String',200);
        set(handles.Nt,'String',800);
        set(handles.Nz,'String',100);
        set(handles.dx,'String',0.25);
        set(handles.dy,'String',0.25);
        set(handles.dz,'String',0.50)
        set(handles.dt_sim,'String',0.05);
        make_grid_Callback(hObject, eventdata, handles);
    case 'HomoCST'
        
        %%%%%%%%%%%%% autocorrelation  %%%%%%%
%         M = evalin('base','M3');
%         template = M(48,:);
%         for i = 1:size(M,1)
%             M2(i,:) = envelope(xcorr(template,M(i,:)));
%             [~, peak(i)] = max(M2(i,:));
%         end
%         Mmin = min(peak);
%         TR = peak*kgrid.dt*1e6;
%         TR = abs(TR-max(TR));
%         assignin('base','TR',TR);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%% Save Images %%%%%%%%%%
        set(handles.savename,'String','LateralFWHM')
        set(handles.fignum,'String','20')
        savefig_Callback(hObject, eventdata, handles)
        set(handles.savename,'String','AxialFWHM')
        set(handles.fignum,'String','19')
        savefig_Callback(hObject, eventdata, handles)
      
       
%         set(handles.ef,'Value',1);
%         set(handles.xrange,'String',num2str([150 250]));
%         set(handles.yrange,'String',num2str([220 320]));
% 
%         pushbutton19_Callback(hObject, eventdata, handles)
%          set(handles.savename,'String','Envelope')
%         set(handles.fignum,'String','9')
%         savefig_Callback(hObject, eventdata, handles)
%         
% %         set(handles.ef,'Value',0);
% %         set(handles.xrange,'String',num2str([1 400]));
% %         set(handles.yrange,'String',num2str([1 400]));
% 
% 
%          set(handles.savename,'String','Lateral')
%         set(handles.fignum,'String','20')
%         savefig_Callback(hObject, eventdata, handles)
%         set(handles.savename,'String','Axial')
%         set(handles.fignum,'String','19')
%         savefig_Callback(hObject, eventdata, handles)
%         
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         set(handles.Nx,'String',256);
%         set(handles.Ny,'String',256);
%         set(handles.Nt,'String',800);
%         set(handles.dx,'String',0.2);
%         set(handles.dy,'String',0.2);
%         set(handles.dt_sim,'String',0.04);
%         make_grid_Callback(hObject, eventdata, handles);
%         set(handles.source_menu1,'Value',4)
%         source_menu1_Callback(hObject, eventdata, handles);
%         set(handles.trans_menu,'Value',6)
%         trans_menu_Callback(hObject, eventdata, handles);
%         set(handles.trans_freq,'String',2)
%         set(handles.trans_cycles,'String',1)
%         set(handles.trans_y0,'String',10)
%         set(handles.trans_focus,'String',num2str([128 170]))
%         set(handles.gauss_apod,'Value',1)
%         set(handles.wave_type,'Value',1)
%         make_medium_Callback(hObject, eventdata, handles);
%         make_source_Callback(hObject, eventdata, handles);
    case 'SingleElements'
        array = str2num(handles.auto_input1.String);
        kgrid = evalin('base','kgrid');
        
        %%%% OLD AUTO FUNCTION
        %         s_array = 102:2:144;
        %         for i = 1:length(array)
        %             set(handles.trans_active_elements,'String',array(i))
        %             make_source_Callback(hObject, eventdata, handles);
        %             set(handles.sim_psim,'Value',0);
        %             simulate_Callback(hObject, eventdata, handles);
        %             for j = 1:length(s_array)
        %                 set(handles.sensormenu,'Value',s_array(j))
        %                 find_max_t_Callback(hObject, eventdata, handles);
        %                 t = str2double(get(handles.tp,'String'));
        %                 T(i,j) = t*(kgrid.dt*1e6);
        %             end
        %         end
        %       assignin('base','T',T);
        
        %%%%% NEW AUTO FUNCTION
        set(handles.show_envelope,'Value',1)
        set(handles.wave_type,'Value',1)
        if handles.auto_cond1.Value
            I_idx = str2num(get(handles.auto_input2,'String'));
            focus = str2num(handles.trans_focus.String);
            sRow = focus(1);
            sCol = focus(2);
            SNUM = kgrid.Nx*(sCol-1)+sRow;
            set(handles.sensormenu,'Value',SNUM)
            for j = 1:length(I_idx)
                new_F = SNUM+I_idx(j);
                for i = 1:length(array)
                    set(handles.trans_active_elements,'String',array(i))
                    make_source_Callback(hObject, eventdata, handles);
                    set(handles.sim_psim,'Value',0);
                    simulate_Callback(hObject, eventdata, handles);
                    set(handles.sensormenu,'Value',new_F)
                    find_max_t_Callback(hObject, eventdata, handles);
                    p = str2double(handles.tp.String);
                    T_P(j,i) = p*(kgrid.dt*1e6);
                    pressure = evalin('base','pressure');
                    pressure = pressure(:,:,p);
                    peak(j,i) = max(pressure(:));
                end
            end
            
        else
            for i = 1:length(array)
                set(handles.trans_active_elements,'String',array(i))
                make_source_Callback(hObject, eventdata, handles);
                set(handles.sim_psim,'Value',0);
                simulate_Callback(hObject, eventdata, handles);
                            pushbutton22_Callback(hObject, eventdata, handles);
                            ae_max_Callback(hObject, eventdata, handles);
                            maxae_Callback(hObject, eventdata, handles);
                            %             pushbutton21_Callback(hObject, eventdata, handles);
                            t = str2double(get(handles.maxtae,'String'));
                            T_AE(i) = t*(kgrid.dt*1e6);
                            M3(i,:) = evalin('base','M');
                %
%                 find_max_t_Callback(hObject, eventdata, handles);
%                 p = str2double(handles.tp.String);
%                 %
%                 T_P(i) = p*(kgrid.dt*1e6);
% 
                pressure = evalin('base','pressure');
                pressure = pressure(:,:,t);
                peak(i) = max(pressure(:));
            end
        end
                assignin('base','T_AE',T_AE);
%         assignin('base','T_P',T_P);
        %         assignin('base','T',T');
        assignin('base','M3',M3);
        sfile = get(handles.auto_input2,'String');
        Apods = peak;
        save(sfile,'M3','Apods','T_AE');
        assignin('base','Apods',peak);
    case 'AETimeReverse'
        array = str2num(handles.auto_input1.String);
        kgrid = evalin('base','kgrid');
        I2 = evalin('base','I_main');
        I_array = [0,round(-256*(1/5.1)),round(256*(1/5.1)),0,...
            round(-256*(2/5.1)),round(256*(2/5.1));0,0,0,round(-256*(1/5.1)),...
            round(-256*(1/5.1)),round(-256*(1/5.1))];
        % s_array = {[128 176],[78 176],[178 176],[128 128],[128 128],[228 128]};
       s_array = str2num(handles.auto_input2.String);
        for i = 1:length(array)
            set(handles.trans_active_elements,'String',array(i))
            make_source_Callback(hObject, eventdata, handles);
            set(handles.sim_psim,'Value',0);
            simulate_Callback(hObject, eventdata, handles);
%             sensor_data = evalin('base','sensor_data');
            pressure = evalin('base','pressure');
%             v = sqrt(size(pressure,1));
%             pressure = reshape(pressure,[v,v,size(pressure,2)]);
%             set(handles.medx,'String',v);
%             set(handles.medy,'String',v);
            assignin('base','pressure',pressure);
            set(handles.use_elec,'Value',1);
            for j = 1:length(I_array)
                focus = s_array{j};
                sRow = focus(1);
                sCol = focus(2);
                SNUM = kgrid.Nx*(sCol-1)+sRow;
                set(handles.sensormenu,'Value',SNUM)
                
                find_max_t_Callback(hObject, eventdata, handles);
                t = str2double(get(handles.tp,'String'));
                T(i,j) = t*(kgrid.dt*1e6);
                I = I2;
                I = circshift(I,I_array(1,j),1);
                I = circshift(I,I_array(2,j),2);
                assignin('base','I',I);
                pushbutton22_Callback(hObject, eventdata, handles);
                ae_max_Callback(hObject, eventdata, handles);
                maxae_Callback(hObject, eventdata, handles);
                pushbutton20_Callback(hObject, eventdata, handles);
                t = str2double(get(handles.maxtae,'String'));
                AE_TR(i,j) = t*(kgrid.dt*1e6);
            end
        end 
        assignin('base','AE_TR',AE_TR);
        assignin('base','TR',T);
    case '3Ddefault'
      
        set(handles.sensor_mode,'Value',4)
          set(handles.Nx,'String',64);
        set(handles.Ny,'String',64);
        set(handles.Nz,'String',64);
        set(handles.Nt,'String',300);
        set(handles.dx,'String',0.5);
        set(handles.dy,'String',0.5);
        set(handles.dz,'String',0.5)
        set(handles.dt_sim,'String',0.05);
        make_grid_Callback(hObject, eventdata, handles);
        set(handles.source_menu1,'Value',4)
        source_menu1_Callback(hObject, eventdata, handles);
        set(handles.trans_menu,'Value',6)
        trans_menu_Callback(hObject, eventdata, handles);
        set(handles.trans_freq,'String',0.5)
        set(handles.trans_cycles,'String',4)
        set(handles.trans_y0,'String',29)
        set(handles.trans_x0,'String',5)
        set(handles.trans_z0,'String',32)
        set(handles.trans_elements,'String',3)
        set(handles.trans_active_elements,'String','1:3')
        set(handles.trans_focus,'String',num2str([64 128]))
        set(handles.gauss_apod,'Value',1)
        set(handles.wave_type,'Value',1)
        set(handles.trans_x_kerf,'String',0.5)
        set(handles.trans_width,'String',1.5)
        set(handles.trans_length,'String',1);
        set(handles.trans_lens_focus,'String',30)
        set(handles.trans_el_focus,'String',45)
        set(handles.use_trans,'Value',1)
        set(handles.sim_3D,'Value',1)
        set(handles.sensor_intensity,'Value',1)
        set(handles.sensor_p_max_all,'Value',0)
        set(handles.sensor_p_max,'Value',0)
        set(handles.sensor_velocity,'Value',0)
        set(handles.sensor_final,'Value',0)
        set(handles.sensor_mode,'Value',4)
        make_sensor_Callback(hObject, eventdata, handles);
        sensor = evalin('base','sensor');
        sensor.mask = ones(64,64,64);
        assignin('base','sensor',sensor);
        
        make_medium_Callback(hObject, eventdata, handles);
        show_tw_Callback(hObject, eventdata, handles);
        make_trans_Callback(hObject,eventdata,handles)
    case 'ConvRecon' %Let TF be a vector of lateral focal points
        TF = str2num(handles.auto_input1.String);
        %         set(handles.wave_type,'Value',2)
        %         set(handles.trans_cycles,'String',3)
        focus = str2num(handles.trans_focus.String);
        fy = focus(2);
        kgrid = evalin('base','kgrid');
        set(handles.pressuremenu,'Value',1)
        for i = 1:length(TF)
            set(handles.trans_focus,'String',num2str([TF(i) focus(2)]))  
            make_source_Callback(hObject, eventdata, handles)
            simulate_Callback(hObject, eventdata, handles);
            pushbutton22_Callback(hObject, eventdata, handles);
            M = evalin('base','M');
            for j = 1:size(M,2)
                M2(:,i,j) = M(:,j);
            end
            
            %             sd = evalin('base','sensor_data');
            %             sd.p = sd.p(:,500:1100);
            %             P_all{i} = gather(reshape(sd.p, [400 400 601]));
            %             source = evalin('base','source');
            %             d1(i) = size(source.p,2);
            %         delays(i) = sqrt((176-10)^2+(128-TF(i))^2);
        end
%         x0 = str2double(handles.trans_x0.String);
%         y0 = str2double(handles.trans_y0.String);
%         for i = 1:length(TF)
%             xdis(i) = TF(i) - x0;
%         end
%         ydis = focus(2) - y0;
%         for i = 1:length(TF)
%             hdis(i) = sqrt(xdis(i)^2+ydis^2);
%         end
%         delays = round(hdis-min(hdis));
%         for i = 1:length(TF)
%             for j = 1:size(M2,3)
%                 M2(:,i,j) = circshift(M2(:,i,j),delays(i)*-1);
%             end
%         end
        %         hdis = hdis*kgrid.dx;
        %         hdist = hdis/1480000;
        %         dels = mod(hdist,kgrid.dt);
        
        
        %         delays = d1-min(d1);
        %         delays = round(delays - min(delays));
        %         assignin('base','delays',delays);
        %         assignin('base','P_all',P_all);
        assignin('base','TF',TF);
        assignin('base','M2',M2);

    case 'HighResTC'
          set(handles.Nx,'String',400);
        set(handles.Ny,'String',400);
        set(handles.Nt,'String',1200);
        set(handles.dx,'String',0.1);
        set(handles.dy,'String',0.1);
        set(handles.dt_sim,'String',0.025);
        set(handles.sensor_mode,'Value',5);
        set(handles.sensor_p_max_all,'Value',0);
        set(handles.sensor_p_max,'Value',0);
        set(handles.sensor_velocity,'Value',0);
        set(handles.sensor_final,'Value',0);
        make_grid_Callback(hObject, eventdata, handles);
        set(handles.source_menu1,'Value',4)
        source_menu1_Callback(hObject, eventdata, handles);
        set(handles.trans_menu,'Value',6)
        trans_menu_Callback(hObject, eventdata, handles);
        set(handles.trans_freq,'String',2)
        set(handles.trans_cycles,'String',3)
        set(handles.trans_y0,'String',10)
        set(handles.trans_focus,'String',num2str([200 270]))
        set(handles.gauss_apod,'Value',1)
        set(handles.wave_type,'Value',2)
        set(handles.tw_encoding,'Value',2);
        make_medium_Callback(hObject, eventdata, handles);
        make_source_Callback(hObject, eventdata, handles);
    case 'Template'
        % Before running, you must get time reversal numbers for each
        % element for a particular medium and then set focus (depth)
        TR = evalin('base','TR');
         %%%%% Include this for apodization %%%%%
         Apods = evalin('base','Apods');
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        source = evalin('base','source');
        
        p = source.p;
        kgrid = evalin('base','kgrid');
        for  i = 1:size(p,1)
            p3(i) = find(p(i,:),1).*kgrid.dt*1e6;
        end
        
        for i = 1:length(TR)
            TR2(i) = abs(TR(i)-max(TR));
            p2(i) = p3(i)-min(p3);
        end
        
        fpts = str2num(get(handles.auto_input1,'String')); %Focal Locations
        elements = get(handles.auto_input2,'String'); %Active Elements
        set(handles.trans_active_elements,'String',elements);
        
        
        template = TR2-p2; %This modification will be applied to all focal pts
        
        %%%% Apply Template to all focal pts %%%%
        focus = str2num(get(handles.trans_focus,'String'));
        z = focus(2);
        x = focus(1);
        if handles.auto_cond2.Value % Runs simulation with template focusng at one pt
            fpt = fpts(1);
            set(handles.trans_focus,'String',[num2str(fpt) ' '  num2str(z)])
            set(handles.usemedia,'Value',0)
            make_source_Callback(hObject, eventdata, handles);
            %%%%% Include this for apodization %%%%%
            set(handles.tr_apod,'Value',1)
            set(handles.tr_all,'Value',1)
            TRbutton_Callback(hObject, eventdata, handles)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            source = evalin('base','source');
            F = source.p;
            for i = 1:size(p,1)
                Pfirst(i) = find(F(i,:),1);
                Plast(i) = find(F(i,:),1,'last');
            end
            Pfirst = min(Pfirst);
            Plast = min(Plast);
            F = padarray(F,[0 1000],0,'both');
            for  i = 1:size(F,1)
                F2(i,:) = circshift(F(i,:),round(template(i)/(kgrid.dt*1e6)));
            end
            for i = 1:size(F2,1)
                Pfirst2(i) = find(F(i,:),1);
                Plast2(i) = find(F(i,:),1,'last');
            end
            Pfirst2 = min(Pfirst2);
            Plast2 = min(Plast2);
            negL = Pfirst2-Pfirst;
            negR = Plast2-Plast+1;
            F2(:,1:negL) = [];
            F2(:,end-negR:end) = [];
            source.p = F2*-1;
            assignin('base','source',source);
            
            simulate_Callback(hObject, eventdata, handles)
        else
            if handles.auto_cond1.Value % Runs analytic analysis of hetero geometric vs homo geometric + template
                
                for i = 1:length(fpts)
                    set(handles.trans_focus,'String',[num2str(fpts(i)) ' ' num2str(z)]);
                    set(handles.usemedia,'Value',1)
                    make_source_Callback(hObject, eventdata, handles);
                    source = evalin('base','source');
                    F = source.p;
                    for  j = 1:size(F,1)
                        F3(j) = find(F(j,:),1)./100; %Heterogenous
                    end
                    for j = 1:length(F3)
                        F2(j) = (F3(j)-min(F3));
                    end
                    
                    
                    set(handles.usemedia,'Value',0)
                    make_source_Callback(hObject, eventdata, handles);
                    source = evalin('base','source');
                    G = source.p;
                    for  j = 1:size(G,1)
                        G3(j) = find(G(j,:),1)./100; %Homogenous
                    end
                    for j = 1:length(G3)
                        G2(j) = G3(j)-(min(G3)+template(j));
                    end
                    
                    Dif(i,:) = F2 - G2;
                    
                end
                assignin('base','Dif',Dif);
            else %Apply the template delays for each focal point and run simulation
                set(handles.pressuremenu,'Value',1)
                for j = 1:length(fpts)
                    set(handles.trans_focus,'String',[num2str(fpts(j)) ' ' num2str(z)]);
                    set(handles.usemedia,'Value',0)
                    make_source_Callback(hObject, eventdata, handles);
                    %%%%% Include this for apodization %%%%%
                    set(handles.tr_apod,'Value',1)
                    set(handles.tr_all,'Value',1)
                    TRbutton_Callback(hObject, eventdata, handles)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    source = evalin('base','source');
                    F = source.p;
                    for i = 1:size(p,1)
                        Pfirst(i) = find(F(i,:),1);
                        Plast(i) = find(F(i,:),1,'last');
                    end
                    Pfirst = min(Pfirst);
                    Plast = min(Plast);
                    F = padarray(F,[0 1000],0,'both');
                    clear F2
                    for  i = 1:size(F,1)
                        F2(i,:) = circshift(F(i,:),round(template(i)/(kgrid.dt*1e6)));
                    end
                    for i = 1:size(F2,1)
                        Pfirst2(i) = find(F(i,:),1);
                        Plast2(i) = find(F(i,:),1,'last');
                    end
                    Pfirst2 = min(Pfirst2);
                    Plast2 = min(Plast2);
                    negL = Pfirst2-Pfirst;
                    negR = Plast2-Plast+1;
                    F2(:,1:negL) = [];
                    F2(:,end-negR:end) = [];
                    source.p = F2*-1;
                    assignin('base','source',source);
                    
                    
                    simulate_Callback(hObject, eventdata, handles)
                    pushbutton22_Callback(hObject, eventdata, handles)
                    M = evalin('base','M');
                    for k = 1:size(M,2)
                        M2(:,j,k) = M(:,k);
                    end
                    
                    %             sd = evalin('base','sensor_data');
                    %             sd.p = sd.p(:,500:1100);
                    %             P_all{i} = gather(reshape(sd.p, [400 400 601]));
                    %             source = evalin('base','source');
                    %             d1(i) = size(source.p,2);
                    %         delays(i) = sqrt((176-10)^2+(128-TF(i))^2);
                end
%                 x0 = str2double(handles.trans_x0.String);
%                 y0 = str2double(handles.trans_y0.String);
%                 for i = 1:length(TR)
%                     xdis(i) = TR(i) - x0;
%                 end
%                 ydis = focus(2) - y0;
%                 for i = 1:length(TR)
%                     hdis(i) = sqrt(xdis(i)^2+ydis^2);
%                 end
%                 delays = round(hdis-min(hdis));
%                 for i = 1:size(M2,2)
%                     for j = 1:size(M2,3)
%                         M2(:,i,j) = circshift(M2(:,i,j),delays(i)*-1);
%                     end
%                 end
                %         hdis = hdis*kgrid.dx;
                %         hdist = hdis/1480000;
                %         dels = mod(hdist,kgrid.dt);
                
                
                %         delays = d1-min(d1);
                %         delays = round(delays - min(delays));
                %         assignin('base','delays',delays);
                %         assignin('base','P_all',P_all);
                assignin('base','TF',TR);
                assignin('base','M2',M2);
            end
        end
        %%%%% FOR SHIFTING AN M2 %%%%%%%
        for i = 1:size(M2,2)
            M2(:,i) = circshift(M2(:,i),-3*round(abs(i-size(M2,2)/2)));
        end
        
        
    case 'TRConv'
        TRfull = evalin('base','TRfull');
        set(handles.wave_type,'Value',1)
        set(handles.tr_all,'Value',1);
        if handles.auto_input1.Value
                APODSfull = evalin('base','APODSfull');
        end
        for i = 1:size(TRfull,1)
            make_source_Callback(hObject, eventdata, handles)
            TR2 = TRfull(i,:);
            assignin('base','TR2',TR2);
            TRbutton_Callback(hObject, eventdata, handles)
            if handles.auto_input1.Value
                set(handles.tr_apod,'Value',1)
                Apods = APODSfull(i,:);
                assignin('base','Apods',Apods);
                TRbutton_Callback(hObject, eventdata, handles)
                set(handles.tr_apod,'Value',0)
            end
            simulate_Callback(hObject, eventdata, handles)
            pushbutton22_Callback(hObject, eventdata, handles)
            M = evalin('base','M');
            for k = 1:size(M,2)
                M2(:,i,k) = M(:,k);
            end
        end
        assignin('base','M2',M2)
    case 'IEEETrans'
        set(handles.Nx,'String',400);
        set(handles.Ny,'String',400);
        set(handles.Nt,'String',1400);
        set(handles.dx,'String',0.1);
        set(handles.dy,'String',0.15);
        set(handles.dt_sim,'String',0.025);
        make_grid_Callback(hObject, eventdata, handles)
        %                 set(handles.trans_menu,'Value',3);
        %                 set(handles.trans_menu,'Value',2);
        set(handles.trans_focus,'String','200 440');
        %                 set(handles.trans_freq,'String',2.5);
        set(handles.sensor_num,'String',117000);
        set(handles.tr_all,'Value',1);
        set(handles.exactx,'String','44');
        set(handles.exacty,'String','0');
        set(handles.exact,'Value',1);
    case '3DfineR'
        set(handles.sensor_mode,'Value',4)
        set(handles.Nx,'String',400);
        set(handles.Ny,'String',400);
        set(handles.Nz,'String',100);
        set(handles.Nt,'String',1200);
        set(handles.dx,'String',0.1);
        set(handles.dy,'String',0.1);
        set(handles.dz,'String',0.2)
        set(handles.dt_sim,'String',0.025);
        make_grid_Callback(hObject, eventdata, handles);
        set(handles.source_menu1,'Value',4)
        source_menu1_Callback(hObject, eventdata, handles);
        set(handles.trans_menu,'Value',2)
        trans_menu_Callback(hObject, eventdata, handles);
        set(handles.trans_freq,'String',2.5)
        set(handles.trans_cycles,'String',3)
        set(handles.trans_y0,'String',200)
        set(handles.trans_x0,'String',20)
        set(handles.trans_z0,'String',50)
        set(handles.trans_elements,'String',96)
        set(handles.az_ele,'String',96)
        set(handles.trans_active_elements,'String','1:96')
        set(handles.trans_focus,'String',num2str([300 200 50]))
        set(handles.gauss_apod,'Value',1)
        set(handles.wave_type,'Value',2)
        set(handles.trans_x_kerf,'String',2)
        set(handles.trans_width,'String',2)
        set(handles.trans_length,'String',1);
        set(handles.trans_lens_focus,'String',30)
        set(handles.trans_el_focus,'String',45)
        set(handles.use_trans,'Value',0)
        set(handles.sim_3D,'Value',1)
        set(handles.use_gpu,'Value',1)
        set(handles.sensor_intensity,'Value',0)
        set(handles.sensor_p_max_all,'Value',0)
        set(handles.sensor_p_max,'Value',0)
        set(handles.sensor_velocity,'Value',0)
        set(handles.sensor_final,'Value',0)
        set(handles.sensor_mode,'Value',4)
        make_sensor_Callback(hObject, eventdata, handles);
        sensor = evalin('base','sensor');
        sensor.record = {'p'};
        sensor.type = {'Custom'};
        set(handles.source_3d,'Value',1)
        
        sensor.mask = ones(64,64,64);
        assignin('base','sensor',sensor);
        make_medium_Callback(hObject, eventdata, handles);
        set(handles.tw_encoding,'Value',2);
        show_tw_Callback(hObject, eventdata, handles);
        make_source_Callback(hObject,eventdata,handles)
    case '3DTest'
        set(handles.sensor_mode,'Value',4)
        set(handles.Nx,'String',100);
        set(handles.Ny,'String',100);
        set(handles.Nz,'String',50);
        set(handles.Nt,'String',500);
        set(handles.dx,'String',1); %Depth
        set(handles.dy,'String',1); %Lateral
        set(handles.dz,'String',1); %Elevational
        set(handles.dt_sim,'String',0.2);
        make_grid_Callback(hObject, eventdata, handles);
        set(handles.source_menu1,'Value',4)
        source_menu1_Callback(hObject, eventdata, handles);
        set(handles.trans_menu,'Value',4)
        trans_menu_Callback(hObject, eventdata, handles);
        set(handles.trans_freq,'String',0.5)
        set(handles.trans_cycles,'String',4)
        set(handles.tw_encoding,'Value',2)
        set(handles.trans_y0,'String',10)
        set(handles.trans_x0,'String',10)
        set(handles.trans_z0,'String',23)
        set(handles.trans_elements,'String',35)
        set(handles.trans_active_elements,'String','1:35')
        set(handles.trans_focus,'String',num2str([50 50]))
        set(handles.gauss_apod,'Value',1)
        set(handles.wave_type,'Value',1)
        set(handles.trans_x_kerf,'String',1)
        set(handles.trans_width,'String',1)
        set(handles.trans_length,'String',4);
        set(handles.trans_lens_focus,'String',30)
        set(handles.trans_el_focus,'String',45)
        set(handles.use_trans,'Value',1)
        set(handles.sim_3D,'Value',1)
        set(handles.sensor_intensity,'Value',0)
        set(handles.sensor_p_max_all,'Value',1)
        set(handles.sensor_p_max,'Value',0)
        set(handles.sensor_velocity,'Value',0)
        set(handles.sensor_final,'Value',0)
        set(handles.sensor_mode,'Value',4)
        make_sensor_Callback(hObject, eventdata, handles);
         set(handles.sensor_mode,'Value',5)
        sensor = evalin('base','sensor');
        sensor.record = {'p'};
        sensor.type = {'Custom'};
        sensor.mask = 1;
        assignin('base','sensor',sensor);
        make_medium_Callback(hObject, eventdata, handles);
        show_tw_Callback(hObject, eventdata, handles);
        make_trans_Callback(hObject,eventdata,handles)
       
                
end

% hObject    handle to autofill (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in automenu.
function automenu_Callback(hObject, eventdata, handles)
H = handles.automenu.String{handles.automenu.Value};
switch H
    case 'HomoCST'
        set(handles.auto_input1,'Visible','off')
        set(handles.auto_cond1,'Visible','off')
        set(handles.auto_cond2,'Visible','off')
        set(handles.auto_input2,'Visible','off')
    case 'SingleElements'
        set(handles.auto_input1,'Visible','on')
        set(handles.auto_cond1,'Visible','on')
        set(handles.auto_cond2,'Visible','off')
        set(handles.auto_input2,'Visible','on')
    case 'AETimeReverse'
         set(handles.auto_input1,'Visible','on')
         set(handles.auto_cond1,'Visible','on')
         set(handles.auto_cond2,'Visible','off')
         set(handles.auto_input2,'Visible','on')
    case 'ConvRecon'
        set(handles.auto_input1,'Visible','on')
        set(handles.auto_cond1,'Visible','off')
        set(handles.auto_cond2,'Visible','off')
        set(handles.auto_input2,'Visible','off')
    case 'Template'
        set(handles.auto_cond1,'Visible','on')
        set(handles.auto_cond2,'Visible','on')
         set(handles.auto_input1,'Visible','on')
         set(handles.auto_input2,'Visible','on')
end
% hObject    handle to automenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns automenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from automenu


% --- Executes during object creation, after setting all properties.
function automenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to automenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ef.
function ef_Callback(hObject, eventdata, handles)
% hObject    handle to ef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ef


% --- Executes on button press in savefig.
function savefig_Callback(hObject, eventdata, handles)
n = str2double(get(handles.fignum,'String'));
figure(n);
d = get(handles.savedir,'String');%'D:\Chet Backup\Manuscripts\Modeling\Neuromod';
name = get(handles.savename,'String');
 set(gca,'Color','none')
export_fig(fullfile(d,name),'-transparent');
% hObject    handle to savefig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function savename_Callback(hObject, eventdata, handles)
% hObject    handle to savename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of savename as text
%        str2double(get(hObject,'String')) returns contents of savename as a double


% --- Executes during object creation, after setting all properties.
function savename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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


% --- Executes on button press in sim_psim.
function sim_psim_Callback(hObject, eventdata, handles)
% hObject    handle to sim_psim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sim_psim



function auto_input1_Callback(hObject, eventdata, handles)
% hObject    handle to auto_input1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of auto_input1 as text
%        str2double(get(hObject,'String')) returns contents of auto_input1 as a double


% --- Executes during object creation, after setting all properties.
function auto_input1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to auto_input1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in maxae.
function maxae_Callback(hObject, eventdata, handles)
V = evalin('base','V');
M = evalin('base','M');
% W = gausswin(55);
if handles.show_envelope.Value
%         for i = 1:size(V,3)
%             V(:,:,i) = abs(hilbert(V(:,:,i)'))';
%     V3 = envelope(M);
    kgrid = evalin('base','kgrid');
%     V3 = ae_demod3(M,kgrid.t_array,0.3*1e6);
%     figure; plot(-V3);
%     title('0.3');
%     V3 = ae_demod3(M,kgrid.t_array,1.25*1e6);
%     figure; plot(V3);
%     title('0.5');
%      V3 = ae_demod3(M,kgrid.t_array,0.4*1e6);
%     figure; plot(V3);
%     title('0.4');
%      V3 = ae_demod3(M,kgrid.t_array,0.35*1e6);
%     figure; plot(V3);
%     title('0.35');
%      V3 = ae_demod3(M,kgrid.t_array,0.45*1e6);
%     figure; plot(V3);
%     title('0.45');
%      V3 = ae_demod3(M,kgrid.t_array,1.3*1e6);
%     figure; plot(V3);
%     title('1.3');
%      V3 = ae_demod3(M,kgrid.t_array,1.5*1e6);
%     figure; plot(V3);
%     title('1.5');
%      V3 = ae_demod3(M,kgrid.t_array,1.7*1e6);
%     figure; plot(V3);
%     title('1.7');
    %         multiWaitbar('Enveloping',i/size(V,3));
%         end
V3 = envelope(M);
    multiWaitbar('CLOSEALL');
else
    V3 = M;
end
% V3 = filtfilt(W,1,V3);
% V2 = max(V,[],1);
% V3 = max(V2,[],2);
T = find(V3 == max(V3),1);
% figure(7); plot(V3);
set(handles.maxtae,'String',T);
% hObject    handle to maxae (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function maxtae_Callback(hObject, eventdata, handles)
% hObject    handle to maxtae (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxtae as text
%        str2double(get(hObject,'String')) returns contents of maxtae as a double


% --- Executes during object creation, after setting all properties.
function maxtae_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxtae (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_movie.
function save_movie_Callback(hObject, eventdata, handles)
% hObject    handle to save_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_movie


% --- Executes on button press in tr_all.
function tr_all_Callback(hObject, eventdata, handles)
% hObject    handle to tr_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tr_all


% --- Executes on button press in resetcolor.
function resetcolor_Callback(hObject, eventdata, handles)
colormap(ones(64,3));
axes(handles.axes1)
hold off
cla reset
axes(handles.axes2)
hold off
cla reset
axes(handles.axes3)
hold off
cla reset
if ismember('L',evalin('base','who'))
    L = evalin('base','L');
    set(handles.LF_Menu,'String',1:size(L,3));
end
if ismember('sensor',evalin('base','who'))
    sensor = evalin('base','sensor');
    set(handles.sensormenu,'String',1:size(sensor.mask,2))
    set(handles.sensormenu,'Visible','on')
    if strcmp(handles.wave_type.String{handles.wave_type.Value},'Focus')
        kgrid = evalin('base','kgrid');
        focus = str2num(handles.trans_focus.String);
        sRow = focus(1);
        sCol = focus(2);
        SNUM = kgrid.Nx/(kgrid.dx*1e4)*(sCol-1)+sRow;
        set(handles.sensormenu,'Value',SNUM)
    end
end
% hObject    handle to resetcolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in plot_metrics.
function plot_metrics_Callback(hObject, eventdata, handles)
% hObject    handle to plot_metrics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_metrics


% --- Executes on button press in select_sensor.
function select_sensor_Callback(hObject, eventdata, handles)
s1 = str2double(get(handles.sensor_1,'String'));
s2 = str2double(get(handles.sensor_2,'String'));
s3 = str2double(get(handles.sensor_3,'String'));
kgrid = evalin('base','kgrid');
sensor = evalin('base','sensor');
pressure = evalin('base','pressure');
if ndims(pressure) == 3
    if handles.box3d.Value
        plane = handles.p_menu.String{handles.p_menu.Value};
        switch plane
            case 'XY'
                q2 = sub2ind([kgrid.Nx kgrid.Ny],s1,s2);
            case 'XZ'
                q2 = sub2ind([kgrid.Nx kgrid.Nz],s1,s2);
            case 'YZ'
                q2 = sub2ind([kgrid.Ny kgrid.Nz],s1,s2);
        end
    else
        q2 = sub2ind([kgrid.Nx kgrid.Ny],s1,s2);
    end
    set(handles.sensormenu,'Value',q2);
else
    q2 = sub2ind([kgrid.Nx kgrid.Ny kgrid.Nz],s1,s2,s3);
    set(handles.sensormenu,'Value',q2);
end

% hObject    handle to select_sensor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function sensor_1_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensor_1 as text
%        str2double(get(hObject,'String')) returns contents of sensor_1 as a double


% --- Executes during object creation, after setting all properties.
function sensor_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensor_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sensor_2_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensor_2 as text
%        str2double(get(hObject,'String')) returns contents of sensor_2 as a double


% --- Executes during object creation, after setting all properties.
function sensor_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensor_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in yline.
function yline_Callback(hObject, eventdata, handles)
% hObject    handle to yline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of yline


% --- Executes on button press in use_gpu.
function use_gpu_Callback(hObject, eventdata, handles)
% hObject    handle to use_gpu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_gpu


% --- Executes on button press in use_trans.
function use_trans_Callback(hObject, eventdata, handles)
% hObject    handle to use_trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_trans


% --- Executes on button press in plot_3d.
function plot_3d_Callback(hObject, eventdata, handles)
% hObject    handle to plot_3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_3d


% --- Executes on button press in Viewer5d.
function Viewer5d_Callback(hObject, eventdata, handles)
beautify;
% hObject    handle to Viewer5d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in sensor_intensity.
function sensor_intensity_Callback(hObject, eventdata, handles)
% hObject    handle to use_conv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in use_conv.
function use_conv_Callback(hObject, eventdata, handles)
% hObject    handle to use_conv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_conv



function savedir_Callback(hObject, eventdata, handles)
% hObject    handle to savedir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of savedir as text
%        str2double(get(hObject,'String')) returns contents of savedir as a double


% --- Executes during object creation, after setting all properties.
function savedir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savedir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function metrics_xr_Callback(hObject, eventdata, handles)
% hObject    handle to metrics_xr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of metrics_xr as text
%        str2double(get(hObject,'String')) returns contents of metrics_xr as a double


% --- Executes during object creation, after setting all properties.
function metrics_xr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to metrics_xr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function metrics_yr_Callback(hObject, eventdata, handles)
% hObject    handle to metrics_yr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of metrics_yr as text
%        str2double(get(hObject,'String')) returns contents of metrics_yr as a double


% --- Executes during object creation, after setting all properties.
function metrics_yr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to metrics_yr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pressuremenu.
function pressuremenu_Callback(hObject, eventdata, handles)
% hObject    handle to pressuremenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pressuremenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pressuremenu


% --- Executes during object creation, after setting all properties.
function pressuremenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pressuremenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in maxP.
function maxP_Callback(hObject, eventdata, handles)
V = abs(evalin('base','pressure'));
V2 = max(V,[],1);
V3 = max(V2,[],2);
V4 = max(V3(200:end));
T = find(V3(100:end) == V4,1)+200;
set(handles.maxtp,'String',T);
t = str2double(handles.tpoint.String);
P = V(:,:,t);
set(handles.dbmax,'String',max(P(:)));
% hObject    handle to maxP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function maxtp_Callback(hObject, eventdata, handles)
% hObject    handle to maxtp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxtp as text
%        str2double(get(hObject,'String')) returns contents of maxtp as a double


% --- Executes during object creation, after setting all properties.
function maxtp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxtp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in exact.
function exact_Callback(hObject, eventdata, handles)
% hObject    handle to exact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of exact


% --- Executes on button press in db_box.
function db_box_Callback(hObject, eventdata, handles)
% hObject    handle to db_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of db_box



function ae_clims_Callback(hObject, eventdata, handles)
% hObject    handle to ae_clims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ae_clims as text
%        str2double(get(hObject,'String')) returns contents of ae_clims as a double


% --- Executes during object creation, after setting all properties.
function ae_clims_CreateFcn(hObject, eventdata, ~)
% hObject    handle to ae_clims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function exactx_Callback(hObject, eventdata, handles)
% hObject    handle to exactx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of exactx as text
%        str2double(get(hObject,'String')) returns contents of exactx as a double


% --- Executes during object creation, after setting all properties.
function exactx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to exactx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function exacty_Callback(hObject, eventdata, handles)
% hObject    handle to exacty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of exacty as text
%        str2double(get(hObject,'String')) returns contents of exacty as a double


% --- Executes during object creation, after setting all properties.
function exacty_CreateFcn(hObject, eventdata, handles)
% hObject    handle to exacty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dbmax_Callback(hObject, eventdata, handles)
% hObject    handle to dbmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dbmax as text
%        str2double(get(hObject,'String')) returns contents of dbmax as a double


% --- Executes during object creation, after setting all properties.
function dbmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dbmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in colorm.
function colorm_Callback(hObject, eventdata, handles)
% hObject    handle to colorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns colorm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from colorm


% --- Executes during object creation, after setting all properties.
function colorm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to colorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in showLF.
function showLF_Callback(hObject, eventdata, handles)
P = evalin('base','L');
L = handles.LF_Menu.Value;

    P = P(:,:,L);

kgrid = evalin('base','kgrid');
xsize = size(P,1);
ysize = size(P,2);
t = size(P,3);
yax = linspace(0,kgrid.dx*xsize,xsize).*1000;
xax = linspace(0,kgrid.dy*ysize,ysize).*1000;
c = size(P,4);
tax = linspace(0,kgrid.dt*t,t)*1000000;
cax = linspace(0,str2num(handles.cdt.String),c);

xrange = str2num(handles.xrange.String);
yrange = str2num(handles.yrange.String);
trange = str2num(handles.trange.String);
crange = str2num(handles.crange.String);

xpt = str2num(handles.xpoint.String);
ypt = str2num(handles.ypoint.String);
tpt = str2num(handles.tpoint.String);
cpt = str2num(handles.cpoint.String);

if handles.ef.Value
    figure(13);
else
axes(handles.axes1);
end
if ~handles.moviebox.Value
    P1 = P(xrange(1):xrange(2),yrange(1):yrange(2));
    imagesc(xax(yrange(1):yrange(2)),yax(xrange(1):xrange(2)),P1,'ButtonDownFcn',{@PlotMetrics,handles});
    caxis([min(P(:)) max(P(:))])
    colormap(ones(64,3));
colormap('hotcold');
    colorbar;
     xlabel('Depth mm')
    ylabel('Lateral mm')
end
% hObject    handle to showLF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in uself.
function uself_Callback(hObject, eventdata, handles)
% hObject    handle to uself (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of uself


% --- Executes on button press in Inversebutton.
function Inversebutton_Callback(hObject, eventdata, handles)
V = evalin('base','V');
P = evalin('base','pressure');
L = evalin('base','L');
kgrid = evalin('base','kgrid');


H = handles.invertmenu.String{handles.invertmenu.Value};

switch H
    case 'demodulation'
        V = permute(V,[1 2 3 5 4]);
        J = zeros(size(V));
        
        for i = 1:size(P,3)
            if handles.uself.Value
                for j = 1:size(L,3)
                    J(:,:,i,j) = (squeeze(V(:,:,i,j)))./L(:,:,j)./P(:,:,i);
                end
            else
                J(:,:,i) = V(:,:,i)./P(:,:,i);
            end
            multiWaitbar('Solving Current Densities',i/size(P,3));
        end
    case 'deconvolution'
        V = squeeze(V(:,:,763,1,:));
        if handles.uself.Value
               dimsL = size(L);
            dimsV = size(V);
              [gridVx gridVy] = meshgrid(linspace(0,1,dimsV(1)),linspace(0,1,dimsV(2)));
            [gridLx gridLy] = meshgrid(linspace(0,1,dimsL(1)),linspace(0,1,dimsL(2)));
            for i = 1:size(V,3)
            V2(:,:,i) = interp2(gridVx,gridVy,V(:,:,i),gridLx,gridLy);
            end
            Vx = zeros(size(V2));
            Vy = Vx;
            Ly = Vx;
            Ly = Vx;
            for k = 1:size(L,3)
                [Vx(:,:,k) Vy(:,:,k)] = gradient(V2(:,:,k));
                [Lx(:,:,k) Ly(:,:,k)] = gradient(L(:,:,k));
            end
           
            
            for i = 1:size(L,1)
                for j = 1:size(L,2)
                    for k = 1:size(L,3)
                        LFM(i,j,k) = Vx(i,j)*Lx(i,j,k)+Vy(i,j)*Ly(i,j,k);
                    end
                end
                multiWaitbar('Generating Lead Field',i/size(L,1));
            end
 
               J = V2./LFM;
            
            
        else %%% NONE OF THIS WORKS
            M2 = evalin('base','M2');
            p = str2double(handles.tpoint.String);
            P = P(180:220,170:370,p);
            for i = 1:size(M2,1)
                for k = 1:size(M2,3)
%                 P2(:,i) = interp1(1:size(P,1),P(:,i),round(1:size(P,1)/2));
%                 P2(:,i) = resample
M3(i,:,k) = resample(M2(i,:,k),2,1);
                end
            end
            for i = 1:size(P,1)
%                 P3(i,:) = interp1(1:size(P,2),P(i,:),round(1:size(P,2)*1.5));
P3(i,:) = resample(P(i,:),3,2);
            end
%             psf = P(:,:,p);
            J = deconvlucy(M3(:,:,1),P3(1:end-1,:)');
        end
        
end
assignin('base','J',J);
multiWaitbar('CLOSEALL');
% hObject    handle to Inversebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in invertmenu.
function invertmenu_Callback(hObject, eventdata, handles)
% hObject    handle to invertmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns invertmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from invertmenu


% --- Executes during object creation, after setting all properties.
function invertmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to invertmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in show_J.
function show_J_Callback(hObject, eventdata, handles)
P = evalin('base','J');
kgrid = evalin('base','kgrid');
xy = size(P,1);
t = size(P,3);
c = size(P,4);
xax = linspace(0,kgrid.dx*xy,xy)*1000;
yax = linspace(0,kgrid.dy*xy,xy)*1000;
tax = linspace(0,kgrid.dt*t,t)*1000000;
cax = linspace(0,str2num(handles.cdt.String),c);

xrange = str2num(handles.xrange.String);
yrange = str2num(handles.yrange.String);
trange = str2num(handles.trange.String);
crange = str2num(handles.crange.String);

xpt = str2num(handles.xpoint.String);
ypt = str2num(handles.ypoint.String);
tpt = str2num(handles.tpoint.String);
cpt = str2num(handles.cpoint.String);

if handles.uself.Value
    L = evalin('base','L');
end
if handles.ef.Value
    figure(11);
else
    axes(handles.axes1);
end
H = handles.invertmenu.String{handles.invertmenu.Value};
if ~handles.moviebox.Value
    switch H
        case 'demodulation'
    if handles.uself.Value
        P1 = P(xrange(1):xrange(2),yrange(1):yrange(2),handles.LF_Menu.Value);
%     P1 = P(xrange(1):xrange(2),yrange(1):yrange(2),tpt(1),handles.LF_Menu.Value);
    else
    P1 = P(xrange(1):xrange(2),yrange(1):yrange(2),tpt(1),cpt(1));
    end
    imagesc(xax(yrange(1):yrange(2)),yax(xrange(1):xrange(2)),P1,'ButtonDownFcn',{@PlotMetrics,handles});
        case 'deconvolution'
            imagesc(P)
    end
    caxis([min(P(:)) max(P(:))])
%     colormap(ones(64,3));
    colormap('hot');
    colorbar;
    xlabel('Depth mm')
    ylabel('Lateral mm')
    title('Current Density');
    
end
if handles.slowtimebox.Value
    tax = cax;
    P2 = squeeze(P(xrange(1):xrange(2),ypt(1),tpt(1),crange(1):crange(2)));
    P3 = squeeze(P(xpt(1),yrange(1):yrange(2),tpt(1),crange(1):crange(2)));
    colorbar;
else
%     P2 = squeeze(P(xrange(1):xrange(2),ypt,trange(1):trange(2),cpt(1)));
%     P3 = squeeze(P(xpt,yrange(1):yrange(2),trange(1):trange(2),cpt(1)));
%     colorbar;
end


% axes(handles.axes2)
% imagesc(tax,xax,P2)
%   xlabel('T us')
%     ylabel('X mm')
% 
% axes(handles.axes3)
% imagesc(tax,yax,P3)
%   xlabel('T us')
%     ylabel('Y mm')

if handles.moviebox.Value
        if handles.save_movie.Value
        clims = [min(P(:)) max(P(:))];
        file = handles.record_filename.String;
        path = handles.record_dir.String;
        v = VideoWriter(fullfile(path,file));
        v.FrameRate = str2double(handles.sim_fps.String);
        open(v);
        p = 'n';
        figure(100)
        T = tpt(1):tpt(2);
        for i = 1:length(T)
            imagesc(P(:,:,T(i)));
            title([p ' = ' num2str(T(i)*kgrid.dt*1e6) ' \mus']);
            colormap('hot');
            text(1,15,[p ' = ' num2str(T(i)*kgrid.dt*1e6) ' \mus'],'Color','white','FontSize',14);  
            caxis(clims)
            xlabel('Depth mm')
            ylabel('Lateral mm')
            frame = getframe;
            writeVideo(v,frame);
        end
        close(v)
    else
    clims = [min(P(:)) max(P(:))];
    axes(handles.axes1);
    if handles.slowtimebox.Value
        T = cpt(1):cpt(2);
        P1 = squeeze(P(xrange(1):xrange(2),yrange(1):yrange(2),tpt(1),crange(1):crange(2)));
    else
    T = tpt(1):tpt(2);
    P1 = squeeze(P(xrange(1):xrange(2),yrange(1):yrange(2),trange(1):trange(2),cpt(1)));
    end
    for i = 1:length(T)
        title(num2str(i));
        imagesc(xax,yax,P1(:,:,T(i)),clims);
        colorbar;
        pause(str2double(handles.movie_pause.String));
    end
    xlabel('Depth mm')
    ylabel('Lateral mm')
        end
end

% hObject    handle to show_J (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in usemedia.
function usemedia_Callback(hObject, eventdata, handles)
% hObject    handle to usemedia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of usemedia


% --- Executes on selection change in LF_Menu.
function LF_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to LF_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns LF_Menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LF_Menu


% --- Executes during object creation, after setting all properties.
function LF_Menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LF_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in auto_cond1.
function auto_cond1_Callback(hObject, eventdata, handles)
% hObject    handle to auto_cond1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of auto_cond1


% --- Executes on button press in auto_cond2.
function auto_cond2_Callback(hObject, eventdata, handles)
% hObject    handle to auto_cond2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of auto_cond2



function auto_input2_Callback(hObject, eventdata, handles)
% hObject    handle to auto_input2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of auto_input2 as text
%        str2double(get(hObject,'String')) returns contents of auto_input2 as a double


% --- Executes during object creation, after setting all properties.
function auto_input2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to auto_input2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tr_apod.
function tr_apod_Callback(hObject, eventdata, handles)
% hObject    handle to tr_apod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tr_apod


% --- Executes on selection change in apodwindow.
function apodwindow_Callback(hObject, eventdata, handles)
% hObject    handle to apodwindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns apodwindow contents as cell array
%        contents{get(hObject,'Value')} returns selected item from apodwindow


% --- Executes during object creation, after setting all properties.
function apodwindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to apodwindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in set_apod.
function set_apod_Callback(hObject, eventdata, handles)
% hObject    handle to set_apod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
H = handles.apodwindow.String{handles.apodwindow.Value};
source = evalin('base','source');
L = size(source.p,2);
L2 = size(source.p2,2);
switch H
    case 'Left'
        source.p2(33:96,:) = zeros(64,L2);
        source.p(33:96,:) = zeros(64,L);
    case 'Right'
        source.p2(1:64,:) = zeros(64,L2);
        source.p(1:64,:) = zeros(64,L);
    case 'Mid'
        source.p2(1:32,:) = zeros(32,L2);
        source.p(1:32,:) = zeros(32,L);
        source.p2(65:96,:) = zeros(32,L2);
        source.p(65:96,:) = zeros(32,L);
    case 'Side'
        source.p2(33:64,:) = zeros(32,L2);
        source.p(33:64,:) = zeros(32,L);
end
assignin('base','source',source)


% --- Executes on button press in shiftapod.
function shiftapod_Callback(hObject, eventdata, handles)
% hObject    handle to shiftapod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
source = evalin('base','source');
S = str2double(handles.apodshift.String);
H = handles.apodwindow.String{handles.apodwindow.Value};
switch H
    case 'Left'
        source.p(1:32,:) = circshift(source.p(1:32,:),S,2);
    case 'Right'
        source.p(65:96,:) = circshift(source.p(65:96,:),S,2);
    case 'Mid'
        source.p(33:64,:) = circshift(source.p(33:64,:),S,2);
    case 'Side'
        source.p(65:96,:) = circshift(source.p(65:96,:),S,2);
         source.p(1:32,:) = circshift(source.p(1:32,:),S,2);
end
assignin('base','source',source);



function apodshift_Callback(hObject, eventdata, handles)
% hObject    handle to apodshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of apodshift as text
%        str2double(get(hObject,'String')) returns contents of apodshift as a double


% --- Executes during object creation, after setting all properties.
function apodshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to apodshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in largesource.
function largesource_Callback(hObject, eventdata, handles)
TR = evalin('base','TR');
kgrid = evalin('base','kgrid');
source = evalin('base','source');
Tx = source.p_mask(:,10);
tx_loc = find(Tx);
I = str2double(handles.currentlength.String);
y = 260;
x1 = tx_loc - 200;
z1 = sqrt(x1.^2+y^2);
for i = 1:length(tx_loc)
    if x1(i) < -I/2
        x2(i) = x1(i)-I/2;
    elseif x1(i) > I/2
         x2(i) = x1(i)+I/2;
    else
        x2(i) = x1(i);
    end
end
z2 = sqrt(x2.^2+y^2);
q = round(z1-z2');
qs = q*kgrid.dt*1e6;
TR2 = TR-qs';
assignin('base','TR2',TR2);
% hObject    handle to largesource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function currentlength_Callback(hObject, eventdata, handles)
% hObject    handle to currentlength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currentlength as text
%        str2double(get(hObject,'String')) returns contents of currentlength as a double


% --- Executes during object creation, after setting all properties.
function currentlength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currentlength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in find_AETR.
function find_AETR_Callback(hObject, eventdata, handles)
M = evalin('base','M3');
% M = envelope(M')';
% if handles.m3_gauss.Value
%     win = gausswin(45);
% end
Q = str2double(handles.M_el.String);
if handles.m3_gauss.Value
    %     tempN = str2double(handles.currentlength.String);
    tempN = Q;
    temp = M(tempN,:);
    M2 = zeros(size(M,1),size(M,2)*2-1);
    for i = 1:size(M,1)
        M2(i,:) = xcorr(temp,M(i,:));
    end
    [~,TR(:)] = max(M2');
     TR = abs(TR-max(TR));
else
    [~,TR(:)] = max(M');
    %     end
end
kgrid = evalin('base','kgrid');
TR = TR.*kgrid.dt.*1e6;
assignin('base','TR2',TR);
% for i = 1:length(TR)-2
%     r = TR(i:i+2);
%     if max(r)-min(r) > 0.4
%         TR(i+2) = median(r);
%     end
% end
% for i = 3:length(TR)-2
%     maxi = length(TR);
%     r = TR(maxi-i-1:maxi-i+1);
%     if max(r)-min(r) > 0.4
%         TR(maxi-i+1) = median(r);
%     end
% end


assignin('base','TR2',TR);

if handles.Mconvert.Value
    M2 = evalin('base','M4');
    n = size(M2,2);
    for i = 1:n
        M = squeeze(M2(:,i,:));
        tempN = Q;
        temp = M(tempN,:);
        M5 = zeros(size(M,1),size(M,2)*2-1);
        for k = 1:size(M,1)
            M5(k,:) = xcorr(temp,M(k,:));
        end
        [~,TR2(:)] = max(M5');
        TR2 = abs(TR2-max(TR2)).*kgrid.dt.*1e6;
        TR3(:,i) = TR2;
    end
    assignin('base','T_AE2',TR3);
end

% source = evalin('base','source');
% I_size = str2double(handles.currentlength.String);
% freq = str2double(handles.trans_freq.String);
% medium = evalin('base','medium');
% c = medium.sound_speed(1,1);
% wav = c/freq/1000; %in mm
% % sdif = I_size/10-wav; 
% sdif = I_size/10;
% mask = source.p_mask;
% masksize = size(mask);
% % col = mod(find(mask,1),masksize(1));
% col = ceil(find(mask,1)/masksize(2));
% us = mask(:,col);
% els = find(us);
% center = masksize(1)/2;
% if mod(masksize(1),2) == 0
%     center = (center+0.5).*kgrid.dy*1000;
% end
% els = els.*kgrid.dy*1e3; %pixel to mm
% num = length(els);
% % for i = 1:num
% %     if abs(center-els(i)) > sdif
% %     end
% % end
% 
% %Try 2
% uslat = abs(els(end)-els(1));
% focus = 26;
% for i = 1:num
%     if abs(center-els(i)) < sdif
%         newx = abs(center-els(i));
%         y = els(i);
%         hNew = sqrt(y^2+newx^2);
%         hOld = focus;
%         hNew = sqrt(focus^2+newx^2);
%         hOld = focus;
%         TR2(i) = TR(i)+(hNew-hOld)/(kgrid.dy*1e3)*(kgrid.dt*1e6);
%     elseif center-els(i) > 0
%         oldx = center-els(i);
%         newx = oldx+sdif;
%         hOld = sqrt(oldx^2+focus^2);
%         hNew = sqrt(newx^2+focus^2);
%         TR2(i) = TR(i)-(hNew-hOld)/(kgrid.dy*1e3)*(kgrid.dt*1e6);
%     elseif center-els(i) < 0
%          oldx = abs(center-els(i));
%         newx = oldx+sdif;
%         hOld = sqrt(oldx^2+focus^2);
%         hNew = sqrt(newx^2+focus^2);
%         TR2(i) = TR(i)-(hNew-hOld)/(kgrid.dy*1e3)*(kgrid.dt*1e6);
%     else
%         TR2(i) = TR(i);
%     end
% end
% TR2 = TR2-min(TR2);
% assignin('base','TR2',TR2);


%         assignin('base','TR2',TR2);




% hObject    handle to find_AETR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in m3_gauss.
function m3_gauss_Callback(hObject, eventdata, handles)
% hObject    handle to m3_gauss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of m3_gauss


% --- Executes on button press in plottr.
function plottr_Callback(hObject, eventdata, handles)
TR = evalin('base','TR2');
TR = TR-min(TR);
figure(24); plot(TR,'k','LineWidth',2.5)
ylabel('delay [/mus]')
xlabel('Element Number')
% hObject    handle to plottr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in run_TR.
function run_TR_Callback(hObject, eventdata, handles)
array = str2num(handles.tr_elements.String);
kgrid = evalin('base','kgrid');

set(handles.show_envelope,'Value',1)
set(handles.wave_type,'Value',1)

Type = handles.tr_list.String{handles.tr_list.Value};
switch Type
    case 'AE'
        for i = 1:length(array)
            set(handles.trans_active_elements,'String',array(i))
            make_source_Callback(hObject, eventdata, handles);
            set(handles.sim_psim,'Value',0);
            simulate_Callback(hObject, eventdata, handles);
            pushbutton22_Callback(hObject, eventdata, handles);
            ae_max_Callback(hObject, eventdata, handles);
            maxae_Callback(hObject, eventdata, handles);
            t = str2double(get(handles.maxtae,'String'));
            T_AE(i) = t*(kgrid.dt*1e6);
            M3(i,:) = evalin('base','M');
            pressure = evalin('base','pressure');
            pressure = pressure(:,:,t);
            peak(i) = max(pressure(:));
        end
        assignin('base','T_AE',T_AE);
        assignin('base','M3',M3);
        Apods = peak;
        if ~isempty(handles.TR_save.String)
            sfile = get(handles.TR_save,'String');
            save(sfile,'M3','Apods','T_AE');
        end
        assignin('base','Apods',peak);
    case 'P'
        if ~isempty(handles.sensor_num.String)
            set(handles.sensormenu,'Value',str2double(handles.sensor_num.String))
        end
        for i = 1:length(array)
            set(handles.trans_active_elements,'String',array(i))
            make_source_Callback(hObject, eventdata, handles);
            set(handles.sim_psim,'Value',0);
            simulate_Callback(hObject, eventdata, handles);
            find_max_t_Callback(hObject, eventdata, handles);
            p = str2double(handles.tp.String);
            T_P(i) = p*(kgrid.dt*1e6);
            pressure = evalin('base','pressure');
            pressure = pressure(:,:,p);
            peak(i) = max(pressure(:));
        end
        assignin('base','T_P',T_P);
        Apods = peak;
        if ~isempty(handles.TR_save.String)
            sfile = get(handles.TR_save,'String');
            save(sfile,'Apods','T_P');
        end
        assignin('base','Apods',peak);
    case 'PR'
        if ~isempty(handles.sensor_num.String)
            set(handles.sensormenu,'Value',str2double(handles.sensor_num.String))
        end
        for i = 1:length(array)
            set(handles.trans_active_elements,'String',array(i))
            make_source_Callback(hObject, eventdata, handles);
            set(handles.sim_psim,'Value',0);
            simulate_Callback(hObject, eventdata, handles);
            find_max_t_Callback(hObject, eventdata, handles);
            p = str2double(handles.tp.String);
            T_P(i) = p*(kgrid.dt*1e6);
            pushbutton22_Callback(hObject, eventdata, handles);
            ae_max_Callback(hObject, eventdata, handles);
            maxae_Callback(hObject, eventdata, handles);
            t = str2double(get(handles.maxtae,'String'));
            T_AE(i) = t*(kgrid.dt*1e6);
            M3(i,:) = evalin('base','M');
            pressure = evalin('base','pressure');
            pressure = pressure(:,:,p);
            peak(i) = max(pressure(:));
        end
        assignin('base','T_P',T_P);
        assignin('base','T_AE',T_AE);
        assignin('base','M3',M3);
        Apods = peak;
        if ~isempty(handles.TR_save.String)
            sfile = get(handles.TR_save,'String');
            save(sfile,'Apods','T_P','M3','T_AE');
        end
        assignin('base','Apods',peak);
        %         errordlg('Not coded yet');
        %         return
    case 'BB'
        if ~isempty(handles.sensor_num.String)
            set(handles.sensormenu,'Value',str2double(handles.sensor_num.String))
        end
        Currents = {'Depth_Dipole_1_0_cm','Depth_Dipole_0_9_cm','Depth_Dipole_0_8_cm',...
            'Depth_Dipole_0_7_cm','Depth_Dipole_0_6_cm','Depth_Dipole_0_5_cm',...
            'Depth_Dipole_0_4_cm','Depth_Dipole_0_3_cm','Depth_Dipole_0_2_cm',...
            'Depth_Dipole_0_1_cm','Depth_Dipole_0_05_cm'};
        for i = 1:length(array)
            set(handles.trans_active_elements,'String',array(i))
            make_source_Callback(hObject, eventdata, handles);
            set(handles.sim_psim,'Value',0);
            simulate_Callback(hObject, eventdata, handles);
            %             find_max_t_Callback(hObject, eventdata, handles);
            %             p = str2double(handles.tp.String);
            %             T_P(i) = p*(kgrid.dt*1e6);
            ae_max_Callback(hObject, eventdata, handles);
            for j = 1:length(Currents)
                load(Currents{j});
                assignin('base','I',I);
                pause(1)
                pushbutton20_Callback(hObject, eventdata, handles);
                pushbutton22_Callback(hObject, eventdata, handles);
                
                %             maxae_Callback(hObject, eventdata, handles);
                M = evalin('base','M');
                Mb = ae_demod3(M,kgrid.t_array,1.5*1e6);
                [~,T_max(i,j)] = max(Mb);
                [~,T_min(i,j)] = min(Mb);
                Mb2 = ae_demod3(M,kgrid.t_array,0.5*1e6);
                [~,T_max2(i,j)] = max(Mb2);
                [~,T_min2(i,j)] = min(Mb2);
                [~,T_max3(i,j)] = max(M);
                [~,T_min3(i,j)] = min(M);
                %                 T_max(i,j) = t*(kgrid.dt*1e6);
                M3(i,j,:) = evalin('base','M');
%                 figure; plot(Mb);
            end
%             pressure = evalin('base','pressure');
%             pressure = pressure(:,:,p);
%             peak(i) = max(pressure(:));
        end
        %          assignin('base','T_P',T_P);
        assignin('base','T_min',T_min);
        assignin('base','T_max',T_max);
        assignin('base','T_min2',T_min2);
        assignin('base','T_max2',T_max2);
        assignin('base','T_min3',T_min3);
        assignin('base','T_max3',T_max3);
        assignin('base','M3',M3);
%         Apods = peak;
        if ~isempty(handles.TR_save.String)
            sfile = get(handles.TR_save,'String');
            save(sfile,'T_max','M3','T_min','T_max2','T_min2','T_min3','T_max3');
        end
%         assignin('base','Apods',peak);
    case 'Monopoles'
         if ~isempty(handles.sensor_num.String)
            set(handles.sensormenu,'Value',str2double(handles.sensor_num.String))
        end
        Currents = {'Monopole_1x1','Monopole_2x2','Monopole_4x4',...
            'Monopole_6x6','Monopole_8x8','Monopole_10x10','Monopole_12x12',...
            'Monopole_14x14','Monopole_16x16','Monopole_18x18','Monopole_20x20'};
        for i = 1:length(array)
            set(handles.trans_active_elements,'String',array(i))
            make_source_Callback(hObject, eventdata, handles);
            set(handles.sim_psim,'Value',0);
            simulate_Callback(hObject, eventdata, handles);
            %             find_max_t_Callback(hObject, eventdata, handles);
            %             p = str2double(handles.tp.String);
            %             T_P(i) = p*(kgrid.dt*1e6);
            ae_max_Callback(hObject, eventdata, handles);
            for j = 1:11
                load(Currents{j});
                assignin('base','I',I);
                pause(1)
                 pushbutton20_Callback(hObject, eventdata, handles);
                pushbutton22_Callback(hObject, eventdata, handles);
                
                %             maxae_Callback(hObject, eventdata, handles);
                M = evalin('base','M');
                Menv = envelope(M);
                Mfilt = filtfilt(gausswin(100),1,envelope(Menv));
                [~,T_AE1(i,j)] = max(Menv);
                [~,T_AE2(i,j)] = max(Mfilt);
%                 [~,T_max(i,j)] = max(Mb);
%                 [~,T_min(i,j)] = min(Mb);
%                     Mb2 = ae_demod3(M,kgrid.t_array,0.5*1e6);
%                 [~,T_max2(i,j)] = max(Mb);
%                 [~,T_min2(i,j)] = min(Mb);
%                 T_max(i,j) = t*(kgrid.dt*1e6);
                M3(i,j,:) = M;
%                 figure; plot(Mb);
            end
%             pressure = evalin('base','pressure');
%             pressure = pressure(:,:,p);
%             peak(i) = max(pressure(:));
        end
        %          assignin('base','T_P',T_P);
        assignin('base','T_AE1',T_AE1);
        assignin('base','T_AE2',T_AE2);
        assignin('base','M3',M3);
%         Apods = peak;
        if ~isempty(handles.TR_save.String)
            sfile = get(handles.TR_save,'String');
            save(sfile,'T_AE1','M3','T_AE2');
        end
end

% hObject    handle to run_TR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function tr_elements_Callback(hObject, eventdata, handles)
% hObject    handle to tr_elements (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tr_elements as text
%        str2double(get(hObject,'String')) returns contents of tr_elements as a double


% --- Executes during object creation, after setting all properties.
function tr_elements_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tr_elements (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in tr_list.
function tr_list_Callback(hObject, eventdata, handles)
% hObject    handle to tr_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tr_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tr_list


% --- Executes during object creation, after setting all properties.
function tr_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tr_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sensor_num_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensor_num as text
%        str2double(get(hObject,'String')) returns contents of sensor_num as a double


% --- Executes during object creation, after setting all properties.
function sensor_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensor_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TR_save_Callback(hObject, eventdata, handles)
% hObject    handle to TR_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TR_save as text
%        str2double(get(hObject,'String')) returns contents of TR_save as a double


% --- Executes during object creation, after setting all properties.
function TR_save_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TR_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function plottr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plottr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in make_TR2.
function make_TR2_Callback(hObject, eventdata, handles)
H = handles.tr2_list.String{handles.tr2_list.Value};
assignin('base','TR2',evalin('base',H));

% hObject    handle to make_TR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in tr2_list.
function tr2_list_Callback(hObject, eventdata, handles)
% hObject    handle to tr2_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tr2_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tr2_list


% --- Executes during object creation, after setting all properties.
function tr2_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tr2_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in comp_box.
function comp_box_Callback(hObject, eventdata, handles)
% hObject    handle to comp_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of comp_box


% --- Executes on button press in template_box.
function template_box_Callback(hObject, eventdata, handles)
% hObject    handle to template_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of template_box


% --- Executes on button press in Mconvert.
function Mconvert_Callback(hObject, eventdata, handles)
% hObject    handle to Mconvert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Mconvert


% --- Executes on button press in Mnums.
function Mnums_Callback(hObject, eventdata, handles)
M4 = evalin('base','M4');
M3 = evalin('base','M3');
M3 = squeeze(M4(:,str2double(handles.M_trial.String),:));
M = M3(str2double(handles.M_el.String),:);
assignin('base','M3',M3);
assignin('base','M',M);
% hObject    handle to Mnums (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in TRnums.
function TRnums_Callback(hObject, eventdata, handles)
T_AE2 = evalin('base','T_AE2');
T_AE = squeeze(T_AE2(:,str2double(handles.M_trial.String),:));
assignin('base','T_AE',T_AE);
% hObject    handle to TRnums (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function M_el_Callback(hObject, eventdata, handles)
% hObject    handle to M_el (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of M_el as text
%        str2double(get(hObject,'String')) returns contents of M_el as a double


% --- Executes during object creation, after setting all properties.
function M_el_CreateFcn(hObject, eventdata, handles)
% hObject    handle to M_el (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function M_trial_Callback(hObject, eventdata, handles)
% hObject    handle to M_trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of M_trial as text
%        str2double(get(hObject,'String')) returns contents of M_trial as a double


% --- Executes during object creation, after setting all properties.
function M_trial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to M_trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in negfocus.
function negfocus_Callback(hObject, eventdata, handles)
% hObject    handle to negfocus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of negfocus


% --- Executes on button press in source_3d.
function source_3d_Callback(hObject, eventdata, handles)
if hObject.Value
    set(handles.text136,'String','Dep Az  El');   
else
    set(handles.text136,'String','Az   Dep');
end
% hObject    handle to source_3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of source_3d



function trans_y_kerf_Callback(hObject, eventdata, handles)
% hObject    handle to trans_y_kerf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_y_kerf as text
%        str2double(get(hObject,'String')) returns contents of trans_y_kerf as a double


% --- Executes during object creation, after setting all properties.
function trans_y_kerf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_y_kerf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function az_ele_Callback(hObject, eventdata, handles)
% hObject    handle to az_ele (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of az_ele as text
%        str2double(get(hObject,'String')) returns contents of az_ele as a double


% --- Executes during object creation, after setting all properties.
function az_ele_CreateFcn(hObject, eventdata, handles)
% hObject    handle to az_ele (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ela_ele_Callback(hObject, eventdata, handles)
% hObject    handle to ela_ele (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ela_ele as text
%        str2double(get(hObject,'String')) returns contents of ela_ele as a double


% --- Executes during object creation, after setting all properties.
function ela_ele_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ela_ele (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sta_Callback(hObject, eventdata, handles)
% hObject    handle to sta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sta as text
%        str2double(get(hObject,'String')) returns contents of sta as a double


% --- Executes during object creation, after setting all properties.
function sta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in hadamard.
function hadamard_Callback(hObject, eventdata, handles)
%set(handles.wave_type,'Value',3)
wtype = handles.hp_menu.String{handles.hp_menu.Value};
switch wtype
    case 'Hadamard'
        Hk = str2double(handles.h_list.String{handles.h_list.Value});
        H_mat = hadamard(str2double(handles.h_list.String{handles.h_list.Value}));
        if handles.auto_cond2.Value
            a = str2double(handles.auto_input2.String);
            atemp = H_mat(:,a);
            atemp2 = H_mat(a,:);
            H_mat(1,:) = atemp;
            H_mat(:,1) = atemp2;
            H_mat(a,:) = ones(1,Hk);

            H_mat(:,a) = ones(1,Hk);
            sum_h2 = sum(H_mat,1);
            sum_h = sum_h2(1);
            switch sum_h
                case 2
                    H_mat(1,1) = -1;
                case -2
                    H_mat(1,1) = 1;
            end
        end
        wave = handles.wave_type.Value;
%         if wave == 2
%             foci = str2num(handles.new_foci.String);
%             f_spots = 1:size(foci,1);
%             c = 1;
%             if mod(Hk,length(f_spots)) ~= 0
%                 errordlg('Number of transmits must be a multiple of number of unique foci');
%             end
%             for i = f_spots
%                 for j = 1:Hk/length(f_spots)
%                     focus{(i-1)*(Hk/length(f_spots))+j} = foci(c,:);
%                 end
%                 c = c+1;
%             end
% 
%         end
        for i = 1:length(H_mat)
            disp(i)
            set(handles.sta,'String',num2str(H_mat(i,:)))
%             if wave == 2
%                 set(handles.trans_focus,'String',num2str(focus{i}))
%                 %                 if i <= Hk/f_spots(1)
%                 %                     set(handles.trans_focus,'String',foci(1,:));
%                 %                 else
%                 %                     set(handles.trans_focus,'String',foci(2,:));
%                 %                 end
%             end
            make_source_Callback(hObject, eventdata, handles);
            simulate_Callback(hObject, eventdata, handles);
            pressure = evalin('base','pressure');
            if ndims(pressure) == 3
                H_pressure(:,:,:,i) = pressure;
            elseif ndims(pressure) == 4
                H_pressure(:,:,:,:,i) = pressure;
            end
        end
    case 'Planes'
        wave = handles.wave_type.Value;
        if wave == 2
            deltas = str2num(handles.new_foci.String);
         %   f_spots = 1:size(foci,1);
            tx = str2double(handles.h_list.String{handles.h_list.Value});
            focus = str2num(handles.trans_focus.String);
            x1 = focus(1);
            z1 = focus(2);
            if deltas(1) == 0
                x_total = zeros(1,tx+1);
            else
                x_total = round((-deltas(1)*tx/2:deltas(1):deltas(1)*tx/2));
            end
            if deltas(2) == 0
                z_total = zeros(1,tx+1);
            else
                z_total = round((-deltas(2)*tx/2:deltas(2):deltas(2)*tx/2));
            end
            x_all = x1+x_total;
            z_all = z1+z_total;
            foci = cat(1,x_all,z_all)';
            for i = 1:tx
                disp(i)
                %  set(handles.sta,'String',num2str(H_mat(i,:)))
                set(handles.trans_focus,'String',num2str(foci(i,:)))
                make_source_Callback(hObject, eventdata, handles);
                simulate_Callback(hObject, eventdata, handles);
                pressure = evalin('base','pressure');
                if ndims(pressure) == 3
                    H_pressure(:,:,:,i) = pressure;
                elseif ndims(pressure) == 4
                    H_pressure(:,:,:,:,i) = pressure;
                end
            end
            set(handles.trans_focus,'String',num2str(focus));
            set(handles.trans_steer,'String',num2str(0));
        else



            %%%%%
            % set(handles.wave_type,'Value',1);
            pdelta = str2double(handles.sta.String);
            pnums = str2double(handles.h_list.String{handles.h_list.Value});
            pcenter = str2double(handles.trans_steer.String);
            ptotal = pdelta*pnums;
            padd = linspace(-ptotal/2,ptotal/2,pnums+1);
            psteer = pcenter+padd;
            for i = 1:pnums+1
                disp(i);
                set(handles.trans_steer,'String',psteer(i))
                make_source_Callback(hObject, eventdata, handles);
                simulate_Callback(hObject, eventdata, handles);
                pressure = evalin('base','pressure');
                if ndims(pressure) == 3
                    H_pressure(:,:,:,i) = pressure;
                elseif ndims(pressure) == 4
                    H_pressure(:,:,:,:,i) = pressure;
                end

            end
        end
    case 'Focus'
        d_foci = str2num(handles.new_foci.String);
        foci = str2num(handles.trans_focus.String);
        n = str2double(handles.h_list.String{handles.h_list.Value})+1;
        if d_foci(1) > 0
            lat_foci = -(n-1)/2*d_foci(1):d_foci(1):(n-1)/2*d_foci(1);
            lat_foci = lat_foci + foci(1);
        else 
            lat_foci(1:n) = ones(n,1).*foci(1);
        end
        if d_foci(2) > 0
            dep_foci = -(n-1)/2*d_foci(2):d_foci(2):(n-1)/2*d_foci(2) + foci(2);
        else 
            dep_foci(1:n) = ones(n,1).*foci(2);
        end
        new_foci = [lat_foci;dep_foci];
        new_foci = new_foci';
        for i = 1:n
            disp(i)
            %  set(handles.sta,'String',num2str(H_mat(i,:)))
            set(handles.trans_focus,'String',num2str(new_foci(i,:)))
            make_source_Callback(hObject, eventdata, handles);
            simulate_Callback(hObject, eventdata, handles);
            pressure = evalin('base','pressure');
            if ndims(pressure) == 3
                H_pressure(:,:,:,i) = pressure;
            elseif ndims(pressure) == 4
                H_pressure(:,:,:,:,i) = pressure;
            end
        end

        set(handles.trans_focus,'String',num2str(foci));

end
if ndims(H_pressure) == 4
    H_pressure = permute(H_pressure,[1 2 4 3]);
elseif ndims(H_pressure) == 5
    H_pressure = permute(H_pressure,[1 2 3 5 4]);
end
assignin('base','H_pressure',H_pressure);
if handles.notify.Value 
    msgbox('All done with Scan!')
end


% hObject    handle to hadamard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in h_list.
function h_list_Callback(hObject, eventdata, handles)
% hObject    handle to h_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns h_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from h_list


% --- Executes during object creation, after setting all properties.
function h_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to h_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in array_2d.
function array_2d_Callback(hObject, eventdata, handles)
% hObject    handle to array_2d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of array_2d



function curve_el_Callback(hObject, eventdata, handles)
% hObject    handle to curve_el (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of curve_el as text
%        str2double(get(hObject,'String')) returns contents of curve_el as a double


% --- Executes during object creation, after setting all properties.
function curve_el_CreateFcn(hObject, eventdata, handles)
% hObject    handle to curve_el (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function noise_var_Callback(hObject, eventdata, handles)
% hObject    handle to noise_var (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noise_var as text
%        str2double(get(hObject,'String')) returns contents of noise_var as a double


% --- Executes during object creation, after setting all properties.
function noise_var_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noise_var (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function noise_snr_Callback(hObject, eventdata, handles)
% hObject    handle to noise_snr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noise_snr as text
%        str2double(get(hObject,'String')) returns contents of noise_snr as a double


% --- Executes during object creation, after setting all properties.
function noise_snr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noise_snr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton55.
function pushbutton55_Callback(hObject, eventdata, handles)
nvar = handles.noise_var.String;
X = evalin('base',nvar);
xmax = max(X(:)).*2; %This *2 is to make up for the -.5 in the next var
noises = (rand(size(X))-.5).*xmax;
snrs = str2double(handles.noise_snr.String);
noise_var = noises./snrs;
clear snrs noises
Xnoise = X + noise_var;
assignin('base',nvar,Xnoise);
% hObject    handle to pushbutton55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in add_noise.
function add_noise_Callback(hObject, eventdata, handles)
% hObject    handle to add_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of add_noise



function had_switch_Callback(hObject, eventdata, handles)
% hObject    handle to had_switch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of had_switch as text
%        str2double(get(hObject,'String')) returns contents of had_switch as a double


% --- Executes during object creation, after setting all properties.
function had_switch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to had_switch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton56.
function pushbutton56_Callback(hObject, eventdata, handles)
P = evalin('base',handles.noise_var.String);
H = hadamard(str2double(handles.h_list.String{handles.h_list.Value}));
n = str2double(handles.had_switch.String);
dims = size(P);
for i = 1:dims(3)
    P(:,:,i,:) = P(:,:,i,:).*H(i,n);
end
assignin('base',get(handles.noise_var,'String'),P)
% hObject    handle to pushbutton56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in had_box.
function had_box_Callback(hObject, eventdata, handles)
if hObject.Value == 0
   %  set(handles.hadamard,'String','Planes')
     set(handles.text138,'String','dTheta (deg)');
     set(handles.text147,'String','Foci')
elseif hObject.Value == 1
    set(handles.text147,'String','Elements per group (3d)')
%    set(handles.hadamard,'String','Hadamard')
     set(handles.text138,'String','STA');
end
% hObject    handle to had_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of had_box



function new_foci_Callback(hObject, eventdata, handles)
% hObject    handle to new_foci (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of new_foci as text
%        str2double(get(hObject,'String')) returns contents of new_foci as a double


% --- Executes during object creation, after setting all properties.
function new_foci_CreateFcn(hObject, eventdata, handles)
% hObject    handle to new_foci (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function el_space_Callback(hObject, eventdata, handles)
% hObject    handle to el_space (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of el_space as text
%        str2double(get(hObject,'String')) returns contents of el_space as a double


% --- Executes during object creation, after setting all properties.
function el_space_CreateFcn(hObject, eventdata, handles)
% hObject    handle to el_space (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function az_space_Callback(hObject, eventdata, handles)
% hObject    handle to az_space (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of az_space as text
%        str2double(get(hObject,'String')) returns contents of az_space as a double


% --- Executes during object creation, after setting all properties.
function az_space_CreateFcn(hObject, eventdata, handles)
% hObject    handle to az_space (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in box3d.
function box3d_Callback(hObject, eventdata, handles)
% hObject    handle to box3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of box3d



function s_plane_Callback(hObject, eventdata, handles)
% hObject    handle to s_plane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of s_plane as text
%        str2double(get(hObject,'String')) returns contents of s_plane as a double


% --- Executes during object creation, after setting all properties.
function s_plane_CreateFcn(hObject, eventdata, handles)
% hObject    handle to s_plane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in p_menu.
function p_menu_Callback(hObject, eventdata, handles)
% hObject    handle to p_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns p_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from p_menu


% --- Executes during object creation, after setting all properties.
function p_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in all_sensor.
function all_sensor_Callback(hObject, eventdata, handles)
if hObject.Value == 0
    set(handles.all_sensor,'String','Plane')
else
    set(handles.all_sensor,'String','All')
end
% hObject    handle to all_sensor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of all_sensor


% --- Executes on button press in reconstruct.
function reconstruct_Callback(hObject, eventdata, handles)
% hObject    handle to reconstruct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of reconstruct



function sensor_3_Callback(hObject, eventdata, handles)
% hObject    handle to sensor_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensor_3 as text
%        str2double(get(hObject,'String')) returns contents of sensor_3 as a double


% --- Executes during object creation, after setting all properties.
function sensor_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensor_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in hp_menu.
function hp_menu_Callback(hObject, eventdata, handles)
w = handles.hp_menu.String{handles.hp_menu.Value};
switch w
    case 'Hadamard'
        set(handles.had_box,'Value',1)
    case {'Planes','Focus'}
        set(handles.had_box,'Value',0)
end
% hObject    handle to hp_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns hp_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from hp_menu


% --- Executes during object creation, after setting all properties.
function hp_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hp_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in notify.
function notify_Callback(hObject, eventdata, handles)
% hObject    handle to notify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of notify
