function varargout = gamut_interactive_app(varargin)
% GAMUT_INTERACTIVE_APP MATLAB code for gamut_interactive_app.fig
%      GAMUT_INTERACTIVE_APP, by itself, creates a new GAMUT_INTERACTIVE_APP or raises the existing
%      singleton*.
%
%      H = GAMUT_INTERACTIVE_APP returns the handle to a new GAMUT_INTERACTIVE_APP or the handle to
%      the existing singleton*.
%
%      GAMUT_INTERACTIVE_APP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GAMUT_INTERACTIVE_APP.M with the given input arguments.
%
%      GAMUT_INTERACTIVE_APP('Property','Value',...) creates a new GAMUT_INTERACTIVE_APP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gamut_interactive_app_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gamut_interactive_app_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gamut_interactive_app

% Last Modified by GUIDE v2.5 30-May-2013 02:38:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gamut_interactive_app_OpeningFcn, ...
                   'gui_OutputFcn',  @gamut_interactive_app_OutputFcn, ...
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


% --- Executes just before gamut_interactive_app is made visible.
function gamut_interactive_app_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gamut_interactive_app (see VARARGIN)

% Need to make a dial for this
handles.use_uplab = true;

handles.rgbgamut = fetch_cielchab_gamut('srgb', 2048, 'face-plus', handles.use_uplab);

if ~isfield(handles.rgbgamut,'lch_chr')
    lch_chr = find_gamut_chr(handles.rgbgamut);
    handles.rgbgamut.lch_chr = lch_chr;
end
max_chr = max(...
    [handles.rgbgamut.lch_chr.Lmin.lch(end,2),...
    handles.rgbgamut.lch_chr.Lmax.lch(end,2)]);
set(handles.chr_slider,'Max',max_chr+1);
set(handles.chr_slider,'Value',50);
set(handles.chr_val,'String',num2str(50));

set(gcf,'Renderer','OpenGL');
set(handles.main_axes,'Color',[0.4663 0.4663 0.4663]);

setup_h_plot(handles);
% set(handles.dim_select,'Value',4);setup_lig_Lab_plot(handles);

% Choose default command line output for gamut_interactive_app
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gamut_interactive_app wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gamut_interactive_app_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% -------------------------------------------------------------------------
% --------------------------- Create functions ---------------------------
% -------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function hue_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hue_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes during object creation, after setting all properties.
function hue_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hue_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function hue_panel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hue_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function lig_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lig_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes during object creation, after setting all properties.
function lig_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lig_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function chr_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chr_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function chr_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chr_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function dim_select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dim_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% Have to manually set this as the GUIDE causes there to be new lines in
% the middle of the names
contents = {'Hue','Lightness (LCHab)','Lightness (Lab)','Chroma','Surface','Mesh', 'Surface (LCh)'};
set(hObject,'String',char(contents));



% -------------------------------------------------------------------------
% --------------------------- Dim selector ---------------------------
% -------------------------------------------------------------------------

% --- Executes on selection change in dim_select.
function dim_select_Callback(hObject, eventdata, handles)
% hObject    handle to dim_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns dim_select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dim_select
contents = cellstr(get(hObject,'String'));
dim_name = contents{get(hObject,'Value')};
switch lower(dim_name)
    case 'hue'
        setup_h_plot(handles)
    case {'lightness','lightness (lab)','lightness (lchab)','lightness (lch)'}
        setup_lig_Lab_plot(handles, dim_name);
    case {'lightness (lch)'}
        setup_lig_Lch_plot(handles);
    case 'chroma'
        setup_c_plot(handles);
    case {'surface','surface (lab)'}
        plot_surface(handles);
    case {'mesh','mesh (lab)'}
        plot_mesh(handles);
    case {'surface lch','surface (lch)'}
        plot_surface_lch(handles);
    otherwise
        error('Unknown dimension setting: %s',dim_name);
end



% -------------------------------------------------------------------------
% --------------------------- Hue plot ---------------------------
% -------------------------------------------------------------------------

function setup_h_plot(handles)

set(handles.hue_panel,'Visible','on');
set(handles.lig_panel,'Visible','off');
set(handles.chr_panel,'Visible','off');

set(handles.main_axes,'NextPlot','replace');
v = str2double(get(handles.hue_val,'String'));
update_h_plot(v,handles);
xlabel('Chroma');
ylabel('Lightness');
grid on;
set(handles.main_axes,'Color',[0.4663 0.4663 0.4663]);
set(handles.main_axes,'XLim',[-150 150],'YLim',[0 100]);
set(handles.main_axes,'NextPlot','replacechildren');


function update_h_plot(v,handles)

v  = round(v);
v  = mod(v,360);
set(handles.hue_slider,'Value',v);
set(handles.hue_val,'String',num2str(v));
v2 = mod(v+180,360);
I1 = handles.rgbgamut.lch(:,3)==v;
I2 = handles.rgbgamut.lch(:,3)==v2;
% plot(...
%     [handles.rgbgamut.lch(I1,2);-flipud(handles.rgbgamut.lch(I2,2))],...
%     [handles.rgbgamut.lch(I1,1); flipud(handles.rgbgamut.lch(I2,1))],...
%     'k');
% set(handles.main_axes,'NextPlot','add');
scatter(...
    [handles.rgbgamut.lch(I1,2);-flipud(handles.rgbgamut.lch(I2,2))],...
    [handles.rgbgamut.lch(I1,1); flipud(handles.rgbgamut.lch(I2,1))],...
    20,...
    [handles.rgbgamut.rgb(I1,:); flipud(handles.rgbgamut.rgb(I2,:))],...
    'fill');
set(handles.main_axes,'NextPlot','replacechildren');


function hue_val_Callback(hObject, eventdata, handles)
% hObject    handle to hue_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hue_val as text
%        str2double(get(hObject,'String')) returns contents of hue_val as a double
v = str2double(get(hObject,'String'));
update_h_plot(v,handles);
guidata(hObject, handles);

% --- Executes on slider movement.
function hue_slider_Callback(hObject, eventdata, handles)
% hObject    handle to hue_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
v = get(hObject,'Value');
update_h_plot(v,handles);
guidata(hObject, handles);



% -------------------------------------------------------------------------
% --------------------------- Lightness plot ---------------------------
% -------------------------------------------------------------------------

function setup_lig_Lch_plot(handles)

set(handles.hue_panel,'Visible','off');
set(handles.lig_panel,'Visible','on');
set(handles.chr_panel,'Visible','off');

set(handles.main_axes,'NextPlot','replace');
v = str2double(get(handles.lig_val,'String'));
update_L_plot(v,handles);

function setup_lig_Lab_plot(handles,dim_name)

set(handles.hue_panel,'Visible','off');
set(handles.lig_panel,'Visible','on');
set(handles.chr_panel,'Visible','off');

cla reset;
set(handles.main_axes,'NextPlot','replace');
% v = str2double(get(handles.lig_val,'String'));
v = get(handles.lig_slider,'Value');
update_L_plot(v,handles);
xlabel('a*');
ylabel('b*');
if strcmp(lower(dim_name),'lightness (lab)')
    grid on;
    set(handles.main_axes,'XTick',-150:25:150)
    set(handles.main_axes,'YTick',-150:25:150)
end
set(handles.main_axes,'Color',[0.4663 0.4663 0.4663]);
set(handles.main_axes,'XLim',[-165 165],'YLim',[-165 165]);
set(handles.main_axes,'NextPlot','replacechildren');

function update_L_plot(v,handles)
v  = round(v*handles.rgbgamut.Lprec)/handles.rgbgamut.Lprec;
set(handles.lig_slider,'Value',v);
set(handles.lig_val,'String',num2str(v));
I1 = handles.rgbgamut.lch(:,1)==v;

contents = cellstr(get(handles.dim_select,'String'));
dim_name = contents{get(handles.dim_select,'Value')};

switch lower(dim_name)
    case 'lightness (lch)'
        polar(handles.main_axes,...
            handles.rgbgamut.lch(I1,3)/360*(2*pi),...
            handles.rgbgamut.lch(I1,2));
        set(handles.main_axes,'Color',[0.4663 0.4663 0.4663]);
    case {'lightness (lab)','lightness','lightness (lchab)'}
% set(handles.main_axes,'NextPlot','replacechildren');
% plot(...
%     handles.rgbgamut.lab(I1,2),...
%     handles.rgbgamut.lab(I1,3),...
%     'k');
% set(handles.main_axes,'NextPlot','add');
        cla;
        set(handles.main_axes,'NextPlot','add');
        if ~strcmp(lower(dim_name),'lightness (lab)')
            add_polar_grid(handles);
        end
        scatter(...
            handles.rgbgamut.lab(I1,2),...
            handles.rgbgamut.lab(I1,3),...
            20,...
            handles.rgbgamut.rgb(I1,:),...
            'fill');
        set(handles.main_axes,'NextPlot','replacechildren');
    otherwise
        dim_name
        error('Unknown dimension setting');
end

function add_polar_grid(handles)
t = 0:pi/32:2*pi;
for r=20:20:120
    plot(r*cos(t),r*sin(t),':','Color',[.2 .2 .2]);
    if mod(r,40)==0
        text(0,r,num2str(r));
    end
end
maxr = 140;
plot(maxr*cos(t),maxr*sin(t),'-','Color',[0 0 0]);
r = [0 maxr];
textr = 155;
for t=0:pi/12:(2*pi-pi/64)
    plot(r*cos(t),r*sin(t),':','Color',[.2 .2 .2]);
    text(textr*cos(t)-7,textr*sin(t),num2str(round(t/(2*pi)*360)));
end


function lig_val_Callback(hObject, eventdata, handles)
% hObject    handle to lig_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lig_val as text
%        str2double(get(hObject,'String')) returns contents of lig_val as a double
v = str2double(get(hObject,'String'));
update_L_plot(v,handles);
guidata(hObject, handles);

% --- Executes on slider movement.
function lig_slider_Callback(hObject, eventdata, handles)
% hObject    handle to lig_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
v = get(hObject,'Value');
update_L_plot(v,handles);
guidata(hObject, handles);



% -------------------------------------------------------------------------
% --------------------------- Chroma plot ---------------------------
% -------------------------------------------------------------------------

function [lch_chr] = find_gamut_chr(g)
max_c = max(g.lch(:,3));
hues = unique(g.lch(:,3));
Lmin.lch = nan((max_c+1)*length(hues),3);
Lmax.lch = nan((max_c+1)*length(hues),3);
i = 0;
for c=0:max_c
    for ihue = 1:length(hues)
        h = hues(ihue);
        sublch = g.lch(g.lch(:,3)==h,:);
        sublch = sublch(sublch(:,2)>=c,:);
        if isempty(sublch); continue; end
        i = i+1;
        [C,I] = min(sublch(:,1));
        Lmin.lch(i,:) = [sublch(I,1) c h];
        [C,I] = max(sublch(:,1));
        Lmax.lch(i,:) = [sublch(I,1) c h];
    end
end
Lmin.lch = Lmin.lch(1:i,:);
Lmax.lch = Lmax.lch(1:i,:);

ga = Lmin.lch(:,2).*cos(Lmin.lch(:,3)/360*(2*pi));
gb = Lmin.lch(:,2).*sin(Lmin.lch(:,3)/360*(2*pi));
Lmin.lab = [Lmin.lch(:,1) ga gb];

ga = Lmax.lch(:,2).*cos(Lmax.lch(:,3)/360*(2*pi));
gb = Lmax.lch(:,2).*sin(Lmax.lch(:,3)/360*(2*pi));
Lmax.lab = [Lmax.lch(:,1) ga gb];

% Move back to RGB so we have a set of colors we can show
Lmin.rgb = gd_lab2rgb(Lmin.lab, handles.use_uplab);
Lmax.rgb = gd_lab2rgb(Lmax.lab, handles.use_uplab);

lch_chr.Lmin = Lmin;
lch_chr.Lmax = Lmax;


function setup_c_plot(handles)

set(handles.hue_panel,'Visible','off');
set(handles.lig_panel,'Visible','off');
set(handles.chr_panel,'Visible','on');

set(handles.main_axes,'NextPlot','replace');
v = str2double(get(handles.chr_val,'String'));
update_c_plot(v,handles);
xlabel('Hue');
ylabel('Lightness');
grid on;
set(handles.main_axes,'Color',[0.4663 0.4663 0.4663]);
set(handles.main_axes,'XLim',[0 360],'YLim',[0 100]);
set(handles.main_axes,'NextPlot','replacechildren');



function update_c_plot(v,handles)

v  = round(v);
set(handles.chr_slider,'Value',v);
set(handles.chr_val,'String',num2str(v));

I = handles.rgbgamut.lch_chr.Lmin.lch(:,2)==v;

% plot(...
%     [handles.rgbgamut.lch(I1,2);-flipud(handles.rgbgamut.lch(I2,2))],...
%     [handles.rgbgamut.lch(I1,1); flipud(handles.rgbgamut.lch(I2,1))],...
%     'k');
% set(handles.main_axes,'NextPlot','add');

scatter(...
    handles.rgbgamut.lch_chr.Lmin.lch(I,3),...
    handles.rgbgamut.lch_chr.Lmin.lch(I,1),...
    20,...
    handles.rgbgamut.lch_chr.Lmin.rgb(I,:),...
    'fill');
set(handles.main_axes,'NextPlot','add');
scatter(...
    handles.rgbgamut.lch_chr.Lmax.lch(I,3),...
    handles.rgbgamut.lch_chr.Lmax.lch(I,1),...
    20,...
    handles.rgbgamut.lch_chr.Lmax.rgb(I,:),...
    'fill');
% plot( 35  , 16.5,'xk')
% plot(110  , 75.9,'xk')
% plot(167.5, 91  ,'xk')
% plot(225  , 75.9,'xk')
% plot(300  , 16.5,'xk')

% plot( [35, 110, 167.5, 225, 300], [16.5, 75.9, 91, 75.4, 16.5], 'ko-');
% plot( [35, 110, 225, 300], [16.5, 90, 75.4, 16.5], 'ko-');

% [xxx,yyy,rgb] = parabola_c_fit(v);
% % [xxx,yyy,rgb] = lin_c_fit(v);
% % plot(xxx,yyy,'-k');
% scatter(...
%     xxx,...
%     yyy,...
%     20,...
%     rgb,...
%     'fill');

set(handles.main_axes,'NextPlot','replacechildren');


function [xxx,yyy,rgb] = parabola_c_fit(c)

% c = 44;
% xx = [35, 167.5, 300];
% yy = [16.5, 90, 16.5];
% c = 44;
% xx = [0 152.5 305];
% yy = [22 90 22];
xx = [0 152.5 305];
yy = [40 82 40];
% c = 36;
% xx = [0, 150, 300];
% yy = [15, 75, 15];

x0 = xx(2); % 167.5;
y0 = yy(2); % 91;

% (y-y0) = m*(x-x0)^2

p4 = (yy(1)-y0) / (xx(1)-x0)^2;

xxx = xx(1):xx(3);
yyy = p4*(xxx-x0).^2 + y0;

L = yyy;
h = xxx;
% c = 44;

a = c.*cos(h/360*(2*pi));
b = c.*sin(h/360*(2*pi));

Lab = [L' a' b'];

rgb = gd_lab2rgb(Lab, handles.use_uplab);


function [xxx,yyy,rgb] = lin_c_fit(c)

% xx = [35, 167.5, 300];
% yy = [16.5, 90, 16.5];
xx = [0, 160, 320];
yy = [15, 80, 15];

xxx = linspace(xx(1),xx(2),100);
xxx = [xxx(1:end-1) linspace(xx(2),xx(3),100)];
yyy = linspace(yy(1),yy(2),100);
yyy = [yyy(1:end-1) linspace(yy(2),yy(3),100)];

L = yyy;
h = xxx;
% c = 44;

a = c.*cos(h/360*(2*pi));
b = c.*sin(h/360*(2*pi));

Lab = [L' a' b'];

rgb = gd_lab2rgb(Lab, handles.use_uplab);


function chr_val_Callback(hObject, eventdata, handles)
% hObject    handle to chr_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of chr_val as text
%        str2double(get(hObject,'String')) returns contents of chr_val as a double
v = str2double(get(hObject,'String'));
update_c_plot(v,handles);
guidata(hObject, handles);


% --- Executes on slider movement.
function chr_slider_Callback(hObject, eventdata, handles)
% hObject    handle to chr_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
v = get(hObject,'Value');
update_c_plot(v,handles);
guidata(hObject, handles);



% -------------------------------------------------------------------------
% --------------------------- Surface plot ---------------------------
% -------------------------------------------------------------------------

function plot_mesh_old(handles)
% Pick resolution which doesn't divide 360 for good coverage
switch handles.rgbgamut.Lprec
    case 1
        res = 7;
    case 2
        res = 13;
    case 3
        res = 23;
    case 4
        res = 29;
    case 5
        res = 37;
    otherwise
%         res = 7*handles.rgbgamut.Lprec;
        res = 5*3*2^(handles.rgbgamut.Lprec-3) -1;
end
cla reset;
set(handles.hue_panel,'Visible','off');
set(handles.lig_panel,'Visible','off');
% fill3(...
%     handles.rgbgamut.lab([1:res:end-1 end],2), ...
%     handles.rgbgamut.lab([1:res:end-1 end],3), ...
%     handles.rgbgamut.lab([1:res:end-1 end],1), ...
%     handles.rgbgamut.lab([1:res:end-1 end],1));
% colormap('');
scatter3(...
    handles.rgbgamut.lab([1:res:end-1 end],2), ...
    handles.rgbgamut.lab([1:res:end-1 end],3), ...
    handles.rgbgamut.lab([1:res:end-1 end],1), ...
    20, ...
    handles.rgbgamut.rgb([1:res:end-1 end],:), ...
    'filled');
set(handles.main_axes,'Color',[0.4663 0.4663 0.4663]);
set(handles.main_axes,'XLim',[-150 150],'YLim',[-150 150],'ZLim',[0 100]);
xlabel('a*')
ylabel('b*')
zlabel('L*')


% -------------------------------------------------------------------------
% --------------------------- Surface plot v2 ---------------------------
% -------------------------------------------------------------------------
function plot_surface(handles)

cla reset;
set(handles.hue_panel,'Visible','off');
set(handles.lig_panel,'Visible','off');

% Get a mesh version of the gamut
if ~isfield(handles.rgbgamut,'lchmesh')
    handles.rgbgamut.lchmesh = make_gamut_mesh(handles.rgbgamut);
end

L = handles.rgbgamut.lchmesh.Lgrid([1:end 1],:);
c = handles.rgbgamut.lchmesh.cgrid([1:end 1],:);
h = handles.rgbgamut.lchmesh.hgrid([1:end 1],:);
a = c.*cosd(h);
b = c.*sind(h);

Lab = [L(:) a(:) b(:)];
CMAP = gd_lab2rgb(Lab, handles.use_uplab);

hs = surf(a,b,L,reshape(CMAP,[size(L) 3]));
set(hs,'EdgeColor','none');
% set(hs,'FaceAlpha',0.75);

set(handles.main_axes,'Color',[0.4663 0.4663 0.4663]);
set(handles.main_axes,'XLim',[-150 150],'YLim',[-150 150],'ZLim',[0 100]);
xlabel('a*')
ylabel('b*')
zlabel('L*')

% -------------------------------------------------------------------------
% --------------------------- Mesh plot -----------------------------------
% -------------------------------------------------------------------------
function plot_mesh(handles)

cla reset;
set(handles.hue_panel,'Visible','off');
set(handles.lig_panel,'Visible','off');

% Get a mesh version of the gamut
if ~isfield(handles.rgbgamut,'lchmesh')
    handles.rgbgamut.lchmesh = make_gamut_mesh(handles.rgbgamut);
end

L = handles.rgbgamut.lchmesh.Lgrid([1:4:(end-1) 1],[1:4:(end-1) end]);
c = handles.rgbgamut.lchmesh.cgrid([1:4:(end-1) 1],[1:4:(end-1) end]);
h = handles.rgbgamut.lchmesh.hgrid([1:4:(end-1) 1],[1:4:(end-1) end]);
a = c.*cosd(h);
b = c.*sind(h);

Lab = [L(:) a(:) b(:)];
CMAP = gd_lab2rgb(Lab, handles.use_uplab);

hs = mesh(a,b,L,reshape(CMAP,[size(L) 3]));
set(hs,'FaceColor','none');

set(handles.main_axes,'Color',[0.4663 0.4663 0.4663]);
set(handles.main_axes,'XLim',[-150 150],'YLim',[-150 150],'ZLim',[0 100]);
xlabel('a*')
ylabel('b*')
zlabel('L*')


% -------------------------------------------------------------------------
% --------------------------- Surface plot v2 ---------------------------
% -------------------------------------------------------------------------
function plot_surface_lch(handles)

cla reset;
set(handles.hue_panel,'Visible','off');
set(handles.lig_panel,'Visible','off');

% Get a mesh version of the gamut
if ~isfield(handles.rgbgamut,'lchmesh')
    handles.rgbgamut.lchmesh = make_gamut_mesh(handles.rgbgamut);
end

L = handles.rgbgamut.lchmesh.Lgrid([1:end 1],:);
c = handles.rgbgamut.lchmesh.cgrid([1:end 1],:);
h = handles.rgbgamut.lchmesh.hgrid([1:end 1],:);
h(end,:) = h(end,:) + 360;
a = c.*cosd(h);
b = c.*sind(h);

Lab = [L(:) a(:) b(:)];
CMAP = gd_lab2rgb(Lab, handles.use_uplab);

hs = surf(h,L,c,reshape(CMAP,[size(L) 3]));
set(hs,'EdgeColor','none');
% set(hs,'FaceAlpha',0.75);

set(handles.main_axes,'Color',[0.4663 0.4663 0.4663]);
set(handles.main_axes,'XLim',[0 360],'YLim',[0 100],'ZLim',[0 150]);
xlabel('h')
ylabel('L')
zlabel('C')


