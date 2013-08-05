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

% Last Modified by GUIDE v2.5 05-Aug-2013 20:36:59

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


handles = fetch_gamut(hObject, eventdata, handles);

% max_chr = max(...
%     [handles.rgbgamut.lch_chr.Lmin.lch(end,2),...
%     handles.rgbgamut.lch_chr.Lmax.lch(end,2)]);
% set(handles.chr_slider,'Max',max_chr+1);
set(handles.chr_slider,'Value',50);
set(handles.chr_val,'String',num2str(50));

set(gcf,'Renderer','OpenGL');
set(handles.main_axes,'Color',[0.4663 0.4663 0.4663]);

hard_refresh_current_plot(handles);

% Choose default command line output for gamut_interactive_app
handles.output = hObject;

% fprintf('%.99f\n',handles.figure1)

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
% contents = {'Hue','Lightness (LCHab)','Lightness (Lab)','Chroma','Surface','Mesh', 'Surface (LCh)'};
contents = {'Hue','Lightness','Chroma','Surface'};
set(hObject,'String',char(contents));


% --- Executes during object creation, after setting all properties.
function plot_intv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_intv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function pcs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pcs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% set(hObject,'SelectionChangeFcn',@displayspace_change);


% -------------------------------------------------------------------------
% --------------------------- Dim selector ---------------------------
% -------------------------------------------------------------------------

% --- Executes on selection change in dim_select.
function dim_select_Callback(hObject, eventdata, handles)
% hObject    handle to dim_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hard_refresh_current_plot(handles)

function hard_refresh_current_plot(handles)

set(handles.hue_panel,'Visible','off');
set(handles.lig_panel,'Visible','off');
set(handles.chr_panel,'Visible','off');
set(handles.surftyp  ,'Visible','off');
set(handles.surfparam,'Visible','off');

% Hints: contents = cellstr(get(hObject,'String')) returns dim_select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dim_select
contents = cellstr(get(handles.dim_select,'String'));
dim_name = contents{get(handles.dim_select,'Value')};
switch lower(dim_name)
    case 'hue'
        setup_h_plot(handles)
        set(handles.hue_panel,'Visible','on');
    case {'lightness','lightness (lab)','lightness (lchab)','lightness (lch)'}
        setup_lig_plot(handles);
        set(handles.lig_panel,'Visible','on');
    case 'chroma'
        setup_c_plot(handles);
        set(handles.chr_panel,'Visible','on');
    case {'surface','surface (lab)','mesh','mesh (lab)','surface lch','surface (lch)'}
        setup_surf_plot(handles);
        set(handles.surftyp  ,'Visible','on');
        set(handles.surfparam,'Visible','on');
        rotate3d on;
    otherwise
        error('Unknown dimension setting: %s',dim_name);
end

function soft_refresh_current_plot(handles)
% Hints: contents = cellstr(get(hObject,'String')) returns dim_select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dim_select
contents = cellstr(get(handles.dim_select,'String'));
dim_name = contents{get(handles.dim_select,'Value')};
switch lower(dim_name)
    case 'hue'
        update_h_plot(handles)
    case {'lightness','lightness (lab)','lightness (lchab)','lightness (lch)'}
        update_L_plot(handles);
    case 'chroma'
        update_c_plot(handles);
    case {'surface','surface (lab)','mesh','mesh (lab)','surface lch','surface (lch)'}
        update_surf(handles);
    otherwise
        error('Unknown dimension setting: %s',dim_name);
end

% -------------------------------------------------------------------------
% --------------------------- Plot interval change  ---------------------------
% -------------------------------------------------------------------------

function plot_intv_Callback(hObject, eventdata, handles)
% hObject    handle to plot_intv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of plot_intv as text
%        str2double(get(hObject,'String')) returns contents of plot_intv as a double

soft_refresh_current_plot(handles);


function intv = get_plot_intv(handles)

intv = str2double(get(handles.plot_intv,'String'));


% -------------------------------------------------------------------------
% --------------------------- Display space change  -----------------------
% -------------------------------------------------------------------------

% --- Executes when selected object is changed in pcs.
function pcs_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in pcs 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

handles = fetch_gamut(hObject, eventdata, handles);
soft_refresh_current_plot(handles);


function handles = fetch_gamut(hObject, eventdata, handles)

switch get(get(handles.pcs,'SelectedObject'),'Tag')
    case 'cielab'
        handles.use_uplab = false;
        handles.cmax      = 150; % 133
    case 'uplab'
        handles.use_uplab = true;
        handles.cmax      = 150; % 150, infinity
    otherwise
        error('Unfamiliar dial position');
end

handles.rgbgamut = fetch_cielchab_gamut('srgb', 2048, 'face-plus', handles.use_uplab);

if ~isfield(handles.rgbgamut,'lch_chr')
    handles.rgbgamut.lch_chr = find_gamut_chr(handles.rgbgamut);
end
if ~isfield(handles.rgbgamut,'lchmesh')
    handles.rgbgamut.lchmesh = make_gamut_mesh(handles.rgbgamut);
end

set(handles.chr_slider,'Max',handles.cmax);

% Update handles structure
guidata(hObject, handles);



% -------------------------------------------------------------------------
% --------------------------- Hue plot ---------------------------
% -------------------------------------------------------------------------

function setup_h_plot(handles)

set(handles.main_axes,'NextPlot','replace');

update_h_plot(handles);

xlabel('Chroma');
ylabel('Lightness');
grid on;
set(handles.main_axes,'Color',[0.4663 0.4663 0.4663]);
set(handles.main_axes,'XLim',[-handles.cmax handles.cmax],'YLim',[0 100]);
set(handles.main_axes,'NextPlot','replacechildren');


function update_h_plot(handles, v)

if nargin<2; v = str2double(get(handles.hue_val,'String')); end

update_h_plot_v2(handles, v);


function update_h_plot_v1(handles, v)

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


function update_h_plot_v2(handles, v)

v  = mod(v,360);
set(handles.hue_slider,'Value',v);
set(handles.hue_val,'String',num2str(v));

% Pixels for every one chroma and one lightness
intv = get_plot_intv(handles);
h = v;
L = 0:intv:100;
c = -handles.cmax:intv:handles.cmax;

cc = meshgrid(c,L);
LL = meshgrid(L,c)';
hh = repmat(h,size(cc));

aa = cc.*cosd(hh);
bb = cc.*sind(hh);
Lab = cat(3,LL,aa,bb);

rgb = gd_lab2rgb(Lab,handles.use_uplab);
li = rgb(:,:,1)<0|rgb(:,:,2)<0|rgb(:,:,3)<0|rgb(:,:,1)>1|rgb(:,:,2)>1|rgb(:,:,3)>1;
rgb(repmat(li,[1 1 3])) = .4663;

image(c,L,rgb);
set(gca,'YDir','normal');

set(handles.main_axes,'NextPlot','replacechildren');


function hue_val_Callback(hObject, eventdata, handles)
% hObject    handle to hue_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hue_val as text
%        str2double(get(hObject,'String')) returns contents of hue_val as a double
v = str2double(get(hObject,'String'));
update_h_plot(handles, v);
guidata(hObject, handles);

% --- Executes on slider movement.
function hue_slider_Callback(hObject, eventdata, handles)
% hObject    handle to hue_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
v = get(hObject,'Value');
v = round(v);
update_h_plot(handles, v);
guidata(hObject, handles);



% -------------------------------------------------------------------------
% --------------------------- Lightness plot ---------------------------
% -------------------------------------------------------------------------

% --- Executes when selected object is changed in lig_param.
function lig_param_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in lig_param 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

setup_lig_plot(handles)


function setup_lig_plot(handles)

cla reset;

switch get(get(handles.lig_param,'SelectedObject'),'Tag')
    case 'lab'
        setup_lig_Lab_plot(handles);
    case 'lch'
        setup_lig_Lch_plot(handles);
    otherwise
        error('Unfamiliar dial position');
end

set(handles.main_axes,'XLim',[-170 170],'YLim',[-170 170]);
set(handles.main_axes,'NextPlot','replacechildren');


function setup_lig_Lch_plot(handles)

update_L_plot(handles);


function setup_lig_Lab_plot(handles)

set(handles.main_axes,'NextPlot','replace');
update_L_plot(handles);
grid on;
set(handles.main_axes,'XTick',-175:25:175);
set(handles.main_axes,'YTick',-175:25:175);
xlabel('a*');
ylabel('b*');
set(handles.main_axes,'Color',[0.4663 0.4663 0.4663]);


function update_L_plot(handles, v)
if nargin<2; v = get(handles.lig_slider,'Value'); end

switch get(get(handles.lig_param,'SelectedObject'),'Tag')
    case 'lab'
        update_L_plot_v2(handles, v);
    case 'lch'
        set(handles.main_axes,'NextPlot','replace');
        update_L_plot_v2_lch(handles, v);
        set(handles.main_axes,'NextPlot','add');
        add_polar_grid(handles);
        set(handles.main_axes,'NextPlot','replacechildren');
        set(handles.main_axes,'XLim',[-170 170],'YLim',[-170 170]);
        axis off;
    otherwise
        error('Unfamiliar dial position');
end


function update_L_plot_v1(handles, v)

v  = round(v/handles.rgbgamut.Lintv)*handles.rgbgamut.Lintv;
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

function update_L_plot_v2(handles, v)

set(handles.lig_slider,'Value',v);
set(handles.lig_val,'String',num2str(v));

intv = get_plot_intv(handles);
L = v;
a = [-fliplr(0:intv:handles.cmax) 0:intv:handles.cmax];
b = a;
aa = meshgrid(a,b);
bb = meshgrid(b,a)';
LL = repmat(L,size(aa));

Lab = cat(3,LL,aa,bb);
rgb = gd_lab2rgb(Lab,handles.use_uplab);
li = rgb(:,:,1)<0|rgb(:,:,2)<0|rgb(:,:,3)<0|rgb(:,:,1)>1|rgb(:,:,2)>1|rgb(:,:,3)>1;
rgb(repmat(li,[1 1 3])) = .4663;

image(a,b,rgb);
set(gca,'YDir','normal');

set(handles.main_axes,'NextPlot','replacechildren');

function update_L_plot_v2_lch(handles, v)

maxc = 160;

set(handles.lig_slider,'Value',v);
set(handles.lig_val,'String',num2str(v));

intv = get_plot_intv(handles);
L = v;
a = [-fliplr(0:intv:maxc) 0:intv:maxc];
b = a;
aa = meshgrid(a,b);
bb = meshgrid(b,a)';
LL = repmat(L,size(aa));

Lab = cat(3,LL,aa,bb);
rgb = gd_lab2rgb(Lab,handles.use_uplab);
li = rgb(:,:,1)<0|rgb(:,:,2)<0|rgb(:,:,3)<0|rgb(:,:,1)>1|rgb(:,:,2)>1|rgb(:,:,3)>1;
rgb(repmat(li,[1 1 3])) = .4663;

li2 = sqrt(aa.^2+bb.^2)>maxc;
bg = get(handles.figure1,'Color');
rgb(repmat(li & li2,[1 1 3])) = mean(bg);

image(a,b,rgb);
set(gca,'YDir','normal');

set(handles.main_axes,'NextPlot','replacechildren');


function add_polar_grid(handles)
t = 0:pi/32:2*pi;
for r=20:20:140
    plot(r*cos(t),r*sin(t),':','Color',[.2 .2 .2]);
    if mod(r,40)==0
        text(0,r,num2str(r));
    end
end
maxr = 160;
plot(maxr*cos(t),maxr*sin(t),'-','Color',[0 0 0]);
r = [0 maxr];
textr = 175;
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
update_L_plot(handles, v);
guidata(hObject, handles);

% --- Executes on slider movement.
function lig_slider_Callback(hObject, eventdata, handles)
% hObject    handle to lig_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
v = get(hObject,'Value');
v = round(v);
update_L_plot(handles, v);
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

set(handles.main_axes,'NextPlot','replace');
update_c_plot(handles);
xlabel('Hue');
ylabel('Lightness');
grid on;
set(handles.main_axes,'Color',[0.4663 0.4663 0.4663]);
set(handles.main_axes,'XLim',[0 360],'YLim',[0 100]);
set(handles.main_axes,'NextPlot','replacechildren');


function update_c_plot(handles, v)
if nargin<2; v = str2double(get(handles.chr_val,'String')); end
update_c_plot_v2(handles,v);


function update_c_plot_v1(handles,v)

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

set(handles.main_axes,'NextPlot','replacechildren');


function update_c_plot_v2(handles, v)

set(handles.chr_slider,'Value',v);
set(handles.chr_val,'String',num2str(v));

intv = get_plot_intv(handles);
c = v;
L = 0:intv:100;
h = 0:intv:360;

hh = meshgrid(h,L);
LL = meshgrid(L,h)';
cc = repmat(c,size(hh));

aa = cc.*cosd(hh);
bb = cc.*sind(hh);
Lab = cat(3,LL,aa,bb);

rgb = gd_lab2rgb(Lab,handles.use_uplab);
li = rgb(:,:,1)<0|rgb(:,:,2)<0|rgb(:,:,3)<0|rgb(:,:,1)>1|rgb(:,:,2)>1|rgb(:,:,3)>1;
rgb(repmat(li,[1 1 3])) = .4663;

image(h,L,rgb);
set(gca,'YDir','normal');

set(handles.main_axes,'NextPlot','replacechildren');

function chr_val_Callback(hObject, eventdata, handles)
% hObject    handle to chr_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of chr_val as text
%        str2double(get(hObject,'String')) returns contents of chr_val as a double
v = str2double(get(hObject,'String'));
update_c_plot(handles, v);
guidata(hObject, handles);


% --- Executes on slider movement.
function chr_slider_Callback(hObject, eventdata, handles)
% hObject    handle to chr_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
v = get(hObject,'Value');
v = round(v);
update_c_plot(handles, v);
guidata(hObject, handles);



% -------------------------------------------------------------------------
% --------------------------- Surface plot ---------------------------
% -------------------------------------------------------------------------

function setup_surf_plot(handles)

cla reset;
update_surf(handles);
view(3);


% --- Executes when selected object is changed in surftyp.
function surftyp_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in surftyp 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

update_surf(handles);
guidata(hObject, handles);

% --- Executes when selected object is changed in surfparam.
function surfparam_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in surfparam 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

update_surf(handles);
guidata(hObject, handles);


function update_surf(handles)

[az,el] = view;
switch get(get(handles.surfparam,'SelectedObject'),'Tag')
    case 'lab'
        switch get(get(handles.surftyp,'SelectedObject'),'Tag')
            case 'surf'
                plot_surface_lab(handles);
            case 'mesh'
                plot_mesh_lab(handles);
            otherwise
                error('Unfamiliar dial position');
        end
    case 'lch'
        switch get(get(handles.surftyp,'SelectedObject'),'Tag')
            case 'surf'
                plot_surface_lch(handles);
            case 'mesh'
                plot_mesh_lch(handles);
            otherwise
                error('Unfamiliar dial position');
        end
    otherwise
        error('Unfamiliar dial position');
end
view(az,el);



% -------------------------------------------------------------------------
% --------------------------- Defunct, scatter mesh -----------------------
% -------------------------------------------------------------------------
function plot_mesh_scatter(handles)
% Pick resolution which doesn't divide 360 for good coverage
switch handles.rgbgamut.Lintv
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
%         res = 7*handles.rgbgamut.Lintv;
        res = 5*3*2^(handles.rgbgamut.Lintv-3) -1;
end
cla reset;
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
set(handles.main_axes,'XLim',[-handles.cmax handles.cmax],'YLim',[-handles.cmax handles.cmax],'ZLim',[0 100]);
xlabel('a*')
ylabel('b*')
zlabel('L*')


% -------------------------------------------------------------------------
% --------------------------- Surface plot v2 ---------------------------
% -------------------------------------------------------------------------
function plot_surface_lab(handles)

L = handles.rgbgamut.lchmesh.Lgrid([1:end 1],:);
c = handles.rgbgamut.lchmesh.cgrid([1:end 1],:);
h = handles.rgbgamut.lchmesh.hgrid([1:end 1],:);
a = c.*cosd(h);
b = c.*sind(h);

Lab = [L(:) a(:) b(:)];
rgb = gd_lab2rgb(Lab, handles.use_uplab);

hs = surf(a,b,L,reshape(rgb,[size(L) 3]));
set(hs,'EdgeColor','none');
% set(hs,'FaceAlpha',0.75);

set(handles.main_axes,'Color',[0.4663 0.4663 0.4663]);
set(handles.main_axes,'XLim',[-handles.cmax handles.cmax],'YLim',[-handles.cmax handles.cmax],'ZLim',[0 100]);
xlabel('a*')
ylabel('b*')
zlabel('L*')

function plot_mesh_lab(handles)

L = handles.rgbgamut.lchmesh.Lgrid([1:4:(end-1) 1],[1:4:(end-1) end]);
c = handles.rgbgamut.lchmesh.cgrid([1:4:(end-1) 1],[1:4:(end-1) end]);
h = handles.rgbgamut.lchmesh.hgrid([1:4:(end-1) 1],[1:4:(end-1) end]);
a = c.*cosd(h);
b = c.*sind(h);

Lab = [L(:) a(:) b(:)];
rgb = gd_lab2rgb(Lab, handles.use_uplab);

hs = mesh(a,b,L,reshape(rgb,[size(L) 3]));
set(hs,'FaceColor','none');

set(handles.main_axes,'Color',[0.4663 0.4663 0.4663]);
set(handles.main_axes,'XLim',[-handles.cmax handles.cmax],'YLim',[-handles.cmax handles.cmax],'ZLim',[0 100]);
xlabel('a*')
ylabel('b*')
zlabel('L*')

function plot_surface_lch(handles)

L = handles.rgbgamut.lchmesh.Lgrid([1:end 1],:);
c = handles.rgbgamut.lchmesh.cgrid([1:end 1],:);
h = handles.rgbgamut.lchmesh.hgrid([1:end 1],:);
h(end,:) = h(end,:) + 360;
a = c.*cosd(h);
b = c.*sind(h);

Lab = [L(:) a(:) b(:)];
rgb = gd_lab2rgb(Lab, handles.use_uplab);

hs = surf(h,L,c,reshape(rgb,[size(L) 3]));
set(hs,'EdgeColor','none');
% set(hs,'FaceAlpha',0.75);

set(handles.main_axes,'Color',[0.4663 0.4663 0.4663]);
set(handles.main_axes,'XLim',[0 360],'YLim',[0 100],'ZLim',[0 handles.cmax]);
xlabel('Hue')
ylabel('Lightness')
zlabel('Chroma')

function plot_mesh_lch(handles)

L = handles.rgbgamut.lchmesh.Lgrid([1:4:(end-1) 1],[1:4:(end-1) end]);
c = handles.rgbgamut.lchmesh.cgrid([1:4:(end-1) 1],[1:4:(end-1) end]);
h = handles.rgbgamut.lchmesh.hgrid([1:4:(end-1) 1],[1:4:(end-1) end]);
h(end,:) = h(end,:) + 360;
a = c.*cosd(h);
b = c.*sind(h);

Lab = [L(:) a(:) b(:)];
rgb = gd_lab2rgb(Lab, handles.use_uplab);

hs = mesh(h,L,c,reshape(rgb,[size(L) 3]));
set(hs,'FaceColor','none');

set(handles.main_axes,'Color',[0.4663 0.4663 0.4663]);
set(handles.main_axes,'XLim',[0 360],'YLim',[0 100],'ZLim',[0 handles.cmax]);
xlabel('Hue')
ylabel('Lightness')
zlabel('Chroma')
