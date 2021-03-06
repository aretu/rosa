function varargout = rosa2019(varargin)
% ROSA2019 M-file for rosa2019.fig
%      ROSA2019, by itself, creates a new ROSA2019 or raises the existing
%      singleton*.
%
%      H = ROSA2019 returns the handle to a new ROSA2019 or the handle to
%      the existing singleton*.
%
%      ROSA2019('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ROSA2019.M with the given input arguments.
%
%      ROSA2019('Property','Value',...) creates a new ROSA2019 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before shivaWIN_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to rosa2019_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help rosa2019

% Last Modified by GUIDE v2.5 23-Oct-2019 19:49:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @rosa2019_OpeningFcn, ...
    'gui_OutputFcn',  @rosa2019_OutputFcn, ...
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

format('long')
% End initialization code - DO NOT EDIT
end

% --- Executes just before rosa2019 is made visible.
function rosa2019_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to rosa2019 (see VARARGIN)
% Function to maximize the window via undocumented Java call.
% Reference: http://undocumentedmatlab.com/blog/minimize-maximize-figure-window
set(handles.figure1, 'units', 'normalized', 'position', [0.01 0.01 0.9 0.9])

% Choose default command line output for rosa2019
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using rosa2019.

% UIWAIT makes rosa2019 wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = rosa2019_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla;

popup_sel_index = get(handles.popupmenu1, 'Value');
switch popup_sel_index
    case 1
        plot(rand(5));
    case 2
        plot(sin(1:0.01:25.99));
    case 3
        bar(1:.5:10);
    case 4
        plot(membrane);
    case 5
        surf(peaks);
end
end


%% WRITE --------------------------------------------------------------------
function write_Callback(hObject, eventdata, handles)
% hObject    handle to write (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[nome,pat]=uiputfile( ...
    {'*.txt', 'All MATLAB Files (*txt)'; ...
    '*.*',                   'All Files (*.*)'}, ...
    'Write as',['~/',handles.filename]);
cd (pat)


h_=findobj('Tag','fluid');
statoF=get(h_,'Value');
h_=findobj('Tag','GH');
statoGH=get(h_,'Value');

if statoF==1
    I={'Time' 'shear1' 'EffPressure' 'Mu1' 'Pf' 'LVDT_low' 'LVDT_high' 'vel' 'slip' 'TempE' 'TempM'};
elseif statoGH==1
    I={'Time' 'shear1' 'Normal' 'Mu1' 'dspring' 'LVDT_low' 'vel' 'slip' 'TempE' 'TempM'};
else
    I={'Time' 'shear1' 'Normal' 'Mu1' 'LVDT_low' 'vel' 'slip' 'TempE'};
end

for j=1:length(I); %1:length(handles.column)
    C(j,1)={'%10.6f '};
    if j==length(I)
        C(j,1)={'%10.6f\n'};
    end
end
C1=cell2mat(C');

for j=1:length(I); %1:length(handles.column)
    N(j,1)={['handles.' I{j} '(l,1),']};
    if j==length(I)
        N(j,1)={['handles.' I{j} '(l,1)']};
    end
end
N1=cell2mat(N');


for j=1:length(I); %1:length(handles.column)
    M(j,1)={['''' I{j} '''' ',']};
    if j==length(I)
        M(j,1)={['''' I{j} '''']};
    end
end
M1=cell2mat(M');

%for j=1:length(handles.column)
%    O(j,1)={['''v' num2str(j) '''' ',']};
%    if j==length(handles.column)
%         O(j,1)={['''v' num2str(j) '''']};
%    end
%end
%O1=cell2mat(O');


for j=1:length(I); %1:length(handles.column)
    S(j,1)={'%s '};
    if j==length(I)
        S(j,1)={'%s\n'};
    end
end
S1=cell2mat(S');

%write in a file
nome2=[nome, 'RED.txt'];
fid = fopen(nome2,'wt');
eval(['fprintf(fid,''' S1 ''',' M1 ');'])

eval(['len=length(handles.' handles.column{1} ');'])

for l=1:len
    eval(['fprintf(fid,''' C1 ''',' N1 ');'])
end
fclose(fid);
if ~ strcmp(fieldnames(handles),'dt'); msgbox(['ATTENTION: handles.dt=none']); end
end

%% OPEN --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%ripulisci precedente

ah_=get(handles.axes1,'children'); set(ah_,'XData',[],'YData',[]);
ah_=get(handles.axes2,'children'); set(ah_,'XData',[],'YData',[]);
ah_=get(handles.axes3,'children'); set(ah_,'XData',[],'YData',[]);

I=strcmp(fieldnames(handles),'column');

if any(I)
    for i=1:length(handles.column)
        eval(['handles=rmfield(handles,''', handles.column{i}, ''');'])
    end
    handles=rmfield(handles,'column');
    
end

fname=fieldnames(handles);
I=strfind(fname,'GEF');
if any(cell2mat(I))
    for i=1:length(fname)
        if ~isempty(strfind(fname{i},'GEF'));
            eval(['handles=rmfield(handles,''', fname{i}, ''');'])
        end
    end
end


handles.load=0; % flag su load o open
clear file1
I=strcmp(fieldnames(handles),'new'); if any(I); handles=rmfield(handles,'new'); end
I=strcmp(fieldnames(handles),'X'); if any(I); handles=rmfield(handles,'X'); end
I=strcmp(fieldnames(handles),'TimeZero'); if any(I); handles=rmfield(handles,'TimeZero'); end



%definisce i grafici da plottare:
%qui ci sono i default
handles.g1=2;
handles.g2=3;
handles.g3=5;

ax_=findobj('Tag','edit1'); set(ax_,'String',handles.g1);
ax_=findobj('Tag','edit2'); set(ax_,'String',handles.g2);
ax_=findobj('Tag','edit3'); set(ax_,'String',handles.g3);

[FileName,PathName] = uigetfile('*.*','All Files (*.*)', ...
    'C:\Users\Stefano\Dropbox\Ricerca\ROSA');


cd (PathName)

%% PART EDITED FOR ROSA
table=readtable(FileName, 'HeaderLines',8);
file1.data = table2array(table(6:end,:));
file1.data=str2double(file1.data);

handles.column=table.Properties.VariableNames;

handles.slope=str2double(table2array(table(3,:))); handles.slope(1,1)=1;
handles.offset=str2double(table2array(table(4,:))); handles.offset(1,1)=0;

handles.filename=FileName;
handles.sm=0;
handles.triggered=0;
handles.cutted=[0 0];
handles.loadT=0;
handles.shearT=0;
ll=1;
nn=size(file1.data,1);

handles.column{1}='Time';
num=length(handles.column);

for n=1:num
    test=double(handles.column{n});
    if any(test==32); handles.column{n}=char(test(test~=32)); end
    handles.(handles.column{n})=file1.data(ll:nn,n);
    %reverse old calibrations?
%     handles.(handles.column{n})=(handles.(handles.column{n})-handles.offset(1,n))/handles.slope(1,n);
    
end

DT=handles.Time(2)-handles.Time(1);

handles.column{num+1}='Stamp';
handles.(handles.column{num+1})= DT*(ones(size(file1.data(ll:nn,1),1),1));

num=length(handles.column);
handles.column{num+1}='Rate';
handles.(handles.column{num+1})= cumsum(ones(size(file1.data(ll:nn,1),1),1));

num=length(handles.column);
handles.column{num+1}='RateZero';
handles.(handles.column{num+1})= cumsum(ones(size(file1.data(ll:nn,1),1),1));

%% --> ele


hv=get(handles.XLab(1),'Value');
handles.TimeZero=cumsum(handles.Stamp);
handles.Time=zeros(size(handles.Stamp));
handles.Time(1)=hv*handles.Stamp(1);
handles.Time(2:end)=hv*handles.Stamp(1) +cumsum(handles.Stamp(2:end)); %plotto il numero di riga

handles.Done=[];

handles.tconv = 1000; %<--- Corrections for time conversion in velocity calculation!

guidata(hObject, handles);
handles.zoom=0;

plotta_ora(handles);
end

function OpenMenuItem2_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%ripulisci precedente

ah_=get(handles.axes1,'children'); set(ah_,'XData',[],'YData',[]);
ah_=get(handles.axes2,'children'); set(ah_,'XData',[],'YData',[]);
ah_=get(handles.axes3,'children'); set(ah_,'XData',[],'YData',[]);

I=strcmp(fieldnames(handles),'column');

if any(I)
    for i=1:length(handles.column)
        eval(['handles=rmfield(handles,''', handles.column{i}, ''');'])
    end
    handles=rmfield(handles,'column');
    
end

fname=fieldnames(handles);
I=strfind(fname,'GEF');
if any(cell2mat(I))
    for i=1:length(fname)
        if ~isempty(strfind(fname{i},'GEF'));
            eval(['handles=rmfield(handles,''', fname{i}, ''');'])
        end
    end
end


handles.load=0; % flag su load o open
clear file1
I=strcmp(fieldnames(handles),'new'); if any(I); handles=rmfield(handles,'new'); end
I=strcmp(fieldnames(handles),'X'); if any(I); handles=rmfield(handles,'X'); end
I=strcmp(fieldnames(handles),'TimeZero'); if any(I); handles=rmfield(handles,'TimeZero'); end



%definisce i grafici da plottare:
%qui ci sono i default
handles.g1=2;
handles.g2=3;
handles.g3=5;

ax_=findobj('Tag','edit1'); set(ax_,'String',handles.g1);
ax_=findobj('Tag','edit2'); set(ax_,'String',handles.g2);
ax_=findobj('Tag','edit3'); set(ax_,'String',handles.g3);

[FileName,PathName] = uigetfile('*.*','All Files (*.*)', ...
    '~/Documents/Roma/Raw lab data');


cd (PathName)

%definisce i parametri da matrice
%handles.column=importdata(FileName,'\t',1);

fid=fopen(FileName,'r');
for i=1:3
    file1=fgets(fid);
end
fclose(fid);


%file0=importdata(FileName,'\t',3);file1=char(file0(3,:));
[I]=find(file1==char(44)); change=false;
if ~isempty(I); file1(I)=char(46); change=true; end
A=sscanf(file1,'%f');
b=length(A); clear A
fid=fopen(FileName,'r');

for i=1:b
    A=fscanf(fid,'%s',1);
    %controlla che non interpreti uno spazio come nuova variabile
    if any(strcmp(fieldnames(handles),'column')) && ...
            strcmp(A,2); handles.column{i-1}={[char(handles.column(i-1)), '2']};
    else
        handles.column(i)={A};
    end
    %controlla che non ce ne siano due uguali
    S=sum(strcmp(handles.column(i), handles.column));
    if S > 1; handles.column{i}=([char(handles.column(i)), '2']); end
    
end

fgets(fid);    fgets(fid); i=0;
if change
    while 1
        i=i+1;
        tline = fgetl(fid);
        if ~ischar(tline); break; end
        [I]=find(tline==char(44));
        if ~isempty(I); tline(I)=char(46); end
        file1.data(i,:)=sscanf(tline,'%f');
    end
else
    file1=importdata(FileName,'\t',3);
end
fclose(fid);

h_=findobj('Tag','dt_value');

%[ndt,vdt]=grp2idx(file1.data(:,1));
%if numel(vdt) > 1; handles.dt=str2double(vdt(2));
%else
%    handles.dt=str2double(vdt(1))
%end


%set(h_,'String',handles.dt);

handles.filename=FileName;

handles.sm=0;
handles.triggered=0;
handles.cutted=[0 0];
handles.loadT=0;
handles.shearT=0;
ll=1;
nn=length(file1.data(:,1));
timess = file1.data(:,1);
if max(timess)>60 || min(timess)>0.7
    disp('Time is in Milliseconds')
    tconv = 1;
elseif max(timess)<60 || min(timess)<0.7
    disp('Time is in seconds')
    tconv = 1000;
else
    disp('Unable to ascertain time units')
end
%primo step:togliere tutto quello che ha un campionamento diverso da dt
%handles.xlab=0:handles.dt:(length(file1.data)-1)*handles.dt;
%memorizza anche gli originali
%eval(['handles.' handles.column{1} ' = cumsum(file1.data(ll:nn,1)); '])
%eval(['handles.v' num2str(1) ' = cumsum(file1.data(ll:nn,1)); '])


handles.column{1}='Time';
num=length(handles.column);

for n=2:num
    test=double(handles.column{n});
    if any(test==32)
        handles.column{n}=char(test(test~=32));
    end
    eval(strcat('handles.',handles.column{n}, '= file1.data(ll:nn,', num2str(n), ');'))
end



handles.column{num+1}='Stamp';
eval(['handles.' handles.column{num+1} '= file1.data(ll:nn,1);'])

num=length(handles.column);
handles.column{num+1}='Rate';
eval(['handles.' handles.column{num+1} '= [1:1:length(file1.data(ll:nn,1))]''; '])

num=length(handles.column);
handles.column{num+1}='RateZero';
eval(['handles.' handles.column{num+1} '= [1:1:length(file1.data(ll:nn,1))]''; '])

% --> ele
hv=get(handles.XLab(1),'Value');
handles.TimeZero=cumsum(handles.Stamp);
handles.Time=zeros(size(handles.Stamp));
handles.Time(1)=hv*handles.Stamp(1);
handles.Time(2:end)=hv*handles.Stamp(1) +cumsum(handles.Stamp(2:end)); %plotto il numero di riga

handles.Done=[];
handles.Time=handles.Time*tconv;
handles.tconv=tconv;

guidata(hObject, handles);
handles.zoom=0;

plotta_ora(handles);
end


%% --- Executes on calling function plotta_ora

function plotta_ora(handles)

hOb=findobj('Tag','XLab');
h_ele=get(hOb,'Value');

if h_ele==2
    handles.X=handles.Time/1000;
elseif h_ele==3;
    if any(strcmp(fieldnames(handles),'slip'))
        handles.X=handles.slip;
    elseif any(strcmp(fieldnames(handles),'Slip_Enc_2'))
        handles.X=handles.Slip_Enc_2;
    end
elseif h_ele==1
    handles.X=handles.Rate;
end



stato=handles.zoom;
posx=get(handles.axes1,'XLim');

eval(['plot(handles.X,handles.' handles.column{(handles.g1)} ',''ob'',''parent'',handles.axes1);']);
legend(handles.axes1,[handles.column{handles.g1}])
lim1=get(handles.axes1,'Ylim'); a=findobj('Tag','lim1S'); set(a,'String',lim1(:,2)); b=findobj('Tag','lim1I'); set(b,'String',lim1(:,1));
if (stato==1); set(handles.axes1,'XLim',[posx]); end

eval(['plot(handles.X,handles.' handles.column{(handles.g2)} ',''ob'',''parent'',handles.axes2);']);
legend(handles.axes2,[handles.column{handles.g2}])
lim2=get(handles.axes2,'Ylim'); a=findobj('Tag','lim2S'); set(a,'String',lim2(:,2)); b=findobj('Tag','lim2I'); set(b,'String',lim2(:,1));
if (stato==1); set(handles.axes2,'XLim',[posx]); end

eval(['plot(handles.X,handles.' handles.column{(handles.g3)} ',''ob'',''parent'',handles.axes3);']);
legend(handles.axes3,[handles.column{handles.g3}])
lim3=get(handles.axes3,'Ylim'); a=findobj('Tag','lim3S'); set(a,'String',lim3(:,2)); b=findobj('Tag','lim3I'); set(b,'String',lim3(:,1));
if (stato ==1) ; set(handles.axes3,'XLim',[posx]); end
end


% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)
end

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
    ['Close ' get(handles.figure1,'Name') '...'],...
    'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});
end


%% --- Executes on button press in zoom.
function zoom_Callback(hObject, eventdata, handles)
% hObject    handle to zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h_ele=get(hObject,'Value');
pippo=[handles.axes1, handles.axes2,handles.axes3];
linkaxes(pippo,'x');


if h_ele==1
    zoom on;
else
    zoom off
end
end
% %zoom xon;
% h_=findobj('Tag','XLab');
% stato=get(h_,'Value');
% if stato==1
% t_cut=handles.Time/1000;
% else
% t_cut=handles.XLab;
% end
%
% [xi,yi]=ginput(2) ;
%
% mat(1,:)=abs(t_cut-ones(size(t_cut))*xi(1));
% mat(2,:)=abs(t_cut-ones(size(t_cut))*xi(2));
% ll(:,1)=find(mat(1,:)==min(mat(1,:)));
% ll(:,2)=find(mat(2,:)==min(mat(2,:)));
%
% %for i=1:length(handles.column)
%
% %eval(['handles.', handles.column{i}, '=handles.' ...
% %    handles.column{i}, '(ll(1,1):ll(1,2),1);'])
% %end
%
% h_=findobj('Tag','XLab');
% stato= get(h_,'Value') ;
% if stato==1; set(h_,'String','time(s)'); handles.X=handles.Time/1000;
% else
% set(h_,'String','xlab'); handles.X= handles.XLab;  %se non voglio anche il xlab triggerato
% end
%
% set(handles.axes1,'Xlim',[handles.X(ll(1,1)) handles.X(ll(1,2))]);
% set(handles.axes2,'Xlim',[handles.X(ll(1,1)) handles.X(ll(1,2))]);
% set(handles.axes3,'Xlim',[handles.X(ll(1,1)) handles.X(ll(1,2))]);
% handles.zoom=1;
%
% guidata(hObject, handles);
%
% plotta_ora(handles)
% linkaxes(pippo,'off');


% --- Executes on button press in pan.
function pan_Callback(hObject, eventdata, handles)
% hObject    handle to pan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h_ele=get(hObject,'Value');
pippo=[handles.axes1, handles.axes2,handles.axes3];
linkaxes(pippo,'x');

if h_ele==1
    pan on;
else
    pan off
end
end

% h_=findobj('Tag','XLab');
% stato= get(h_,'Value') ;
% if stato==1
%     set(h_,'String','time(s)')
% handles.X=handles.Time/1000;
% else
% set(h_,'String','xlab')
% handles.X= handles.XLab;  %se non voglio anche il xlab triggerato
% end
% set(handles.axes1,'XLim',[handles.X(1) handles.X(end)]);
% set(handles.axes2,'XLim',[handles.X(1) handles.X(end)]);
% set(handles.axes3,'XLim',[handles.X(1) handles.X(end)]);
%
% handles.zoom=0;
% plotta_ora(handles)
% guidata(hObject, handles);


%% EDIT1
function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.g1=str2double(get(hObject,'String'));

if handles.g1 ~= handles.g1
    [s,v] = listdlg('PromptString','Select a file:',...
        'SelectionMode','single',...
        'ListString',handles.column);
    
    handles.g1=s; %str2double(get(hObject,'String'));
end

guidata(hObject, handles);
set(hObject,'String',handles.g1)

hOb=findobj('Tag','XLab');
h_ele=get(hOb,'Value');

if h_ele==2
    handles.X=handles.Time/1000;
elseif h_ele==3;
    if any(strcmp(fieldnames(handles),'slip'))
        handles.X=handles.slip;
    elseif any(strcmp(fieldnames(handles),'Slip_Enc_2'))
        handles.X=handles.Slip_Enc_2;
    end
elseif h_ele==1
    handles.X=handles.Rate;
end


eval(['plot(handles.X,handles.' handles.column{(handles.g1)} ',''ob'',''parent'',handles.axes1);']);
legend(handles.axes1,[handles.column{handles.g1}])
end

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
end


%% EDIT2
function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
handles.g2=str2double(get(hObject,'String'));

if handles.g2 ~= handles.g2
    [s,v] = listdlg('PromptString','Select a file:',...
        'SelectionMode','single',...
        'ListString',handles.column);
    
    handles.g2=s; %str2double(get(hObject,'String'));
end
guidata(hObject, handles);
set(hObject,'String',handles.g2)

hOb=findobj('Tag','XLab');
h_ele=get(hOb,'Value');

if h_ele==2
    handles.X=handles.Time/1000;
elseif h_ele==3;
    if any(strcmp(fieldnames(handles),'slip'))
        handles.X=handles.slip;
    elseif any(strcmp(fieldnames(handles),'Slip_Enc_2'))
        handles.X=handles.Slip_Enc_2;
    end
elseif h_ele==1
    handles.X=handles.Rate;
end


eval(['plot(handles.X,handles.' handles.column{(handles.g2)} ',''ob'',''parent'',handles.axes2);']);
legend(handles.axes2,[handles.column{handles.g2}])
end

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


%% EDIT3
function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
%handles.g3=str2double(get(hObject,'String'));
handles.g3=str2double(get(hObject,'String'));

if handles.g3 ~= handles.g3
    [s,v] = listdlg('PromptString','Select a file:',...
        'SelectionMode','single',...
        'ListString',handles.column);
    
    handles.g3=s; %str2double(get(hObject,'String'));
end
guidata(hObject, handles);
set(hObject,'String',handles.g3)


hOb=findobj('Tag','XLab');
h_ele=get(hOb,'Value');

if h_ele==2
    handles.X=handles.Time/1000;
elseif h_ele==3;
    if any(strcmp(fieldnames(handles),'slip'))
        handles.X=handles.slip;
    elseif any(strcmp(fieldnames(handles),'Slip_Enc_2'))
        handles.X=handles.Slip_Enc_2;
    end
elseif h_ele==1
    handles.X=handles.Rate;
end
eval(['plot(handles.X,handles.' handles.column{(handles.g3)} ',''ob'',''parent'',handles.axes3);']);
legend(handles.axes3,[handles.column{handles.g3}])
end

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end




% --- Executes on button press in XLab.
function XLab_Callback(hObject, eventdata, handles)
% hObject    handle to XLab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% hObject    handle to XLab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of XLab

h_ele=get(hObject,'Value')

if h_ele==2
    handles.X=handles.Time/1000;
elseif h_ele==3;
    if any(strcmp(fieldnames(handles),'SlipPC'))
        handles.X=handles.SlipPC;
        %     elseif any(strcmp(fieldnames(handles),'Slip_Enc_2'))
        %         handles.X=handles.Slip_Enc_2;
    else disp('not existent field Slip'); set(hObject,'Value',1);
    end
elseif h_ele==4;
    if any(strcmp(fieldnames(handles),'SlipAIncr'))
        handles.X=handles.SlipAIncr;
        %     elseif any(strcmp(fieldnames(handles),'Slip_Enc_2'))
        %         handles.X=handles.Slip_Enc_2;
    else disp('not existent field Slip'); set(hObject,'Value',1);
    end
elseif h_ele==1
    handles.X=handles.Rate;
end


guidata(hObject, handles);

plotta_ora(handles)
end


% --- Executes on button press in offset.
function offset_Callback(hObject, eventdata, handles)
% hObject    handle to offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to offset_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[sel,v] = listdlg('PromptString','Select a file:',...
    'SelectionMode','multiple',...
    'ListString',handles.column);

h_=handles.column(sel);
hOb=findobj('Tag','XLab');
h_ele=get(hOb,'Value');
htype=get(hOb,'String');

if  strcmp(htype{h_ele},'Time (s)')
    t_cut=handles.Time/1000;
elseif  strcmp(htype{h_ele},'Slip (m)')
    if any(strcmp(fieldnames(handles),'slip'))
        t_cut=handles.slip;
    elseif any(strcmp(fieldnames(handles),'Slip_Enc_2'))
        t_cut=handles.Slip_Enc_2;
    end
else
    
    t_cut=handles.X;
end



[xi,yi]=ginput(1) ;
mat(1,:)=abs(t_cut-ones(size(t_cut))*xi(1));
ll(:,1)=find(mat(1,:)==min(mat(1,:)));

ax_=get(gcf,'CurrentAxes');
for n=1:3
    eval(['s=find(ax_==handles.axes' num2str(n) ');'])
    if s==1; break; end
end
posy=get(ax_,'Ylim');
posx=get(ax_,'Xlim');

eval(['h_=findobj(''Tag'',''edit' num2str(n) ''');']); %n=numero asse
colonna=str2double(get(h_,'String'));                        %numero della colonna

for n=sel
    eval(['handles.' handles.column{n} ' =handles.' handles.column{n} '-handles.' handles.column{n} '(ll(:,1),:);'])
    %eval(['handles.' handles.column{n} '(1:ll(:,1),:)=0;'])
end

nsel=find(strcmp(handles.column(sel),'Axial'));
if isempty(nsel)
    handles.shearT=handles.RateZero(ll);
else
    handles.loadT=handles.RateZero(ll);
end

h_=findobj('Tag','edit1');
s1=str2double(get(h_,'String'));
h_=findobj('Tag','edit2');
s2=str2double(get(h_,'String'));
h_=findobj('Tag','edit3');
s3=str2double(get(h_,'String'));

guidata(hObject, handles);
plotta_ora(handles)
end




% --- Executes on button press in trigger.
function trigger_Callback(hObject, eventdata, handles)
% hObject    handle to trigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hOb=findobj('Tag','XLab');
h_ele=get(hOb,'Value');

if h_ele==2
    t_cut=handles.Time/1000;
elseif h_ele==3;
    if any(strcmp(fieldnames(handles),'slip'))
        t_cut=handles.slip;
    elseif any(strcmp(fieldnames(handles),'Slip_Enc_2'))
        t_cut=handles.Slip_Enc_2;
    end
elseif h_ele==1
    t_cut=handles.XLab;
end


A=find(strcmp(fieldnames(handles),'shearT'));
if isempty(A)
    [xi,yi]=ginput(1) ;
    mat(1,:)=abs(t_cut-ones(size(t_cut))*xi(1));
    ll(:,1)=find(mat(1,:)==min(mat(1,:)));
else
    prev_trig=find(handles.RateZero==handles.shearT);
    k=menu(['trigger is' num2str(handles.shearT) '. Is that ok?'],'si','no','man');
    
    if (k==1)
        ll(:,1)=prev_trig;
    elseif (k==2)
        [xi,yi]=ginput(1) ;
        mat(1,:)=abs(t_cut-ones(size(t_cut))*xi(1));
        ll(:,1)=find(mat(1,:)==min(mat(1,:)));
    elseif (k==3)
        t_cut=handles.XLab;
        [xi]=input('digit a triggering number here = ') ;
        mat(1,:)=abs(t_cut-ones(size(t_cut))*xi(1));
        ll(:,1)=find(mat(1,:)==min(mat(1,:)));
        
    end
end
handles.triggered=handles.RateZero(ll);

%cycle if running in t=0
%for n=1:length(handles.column)
I(1)=find(strcmp(handles.column,'Rate'));

%for n=I
%eval(['handles.' handles.column{n} ' =handles.' handles.column{n} '-handles.' handles.column{n} '(ll(:,1),:);'])
%endhandle
handles.Time=handles.Time-handles.Time(ll(:,1),:);

set(hOb,'Value',2)


ax_=get(handles.axes1,'Children'); dataY=get(ax_,'YData'); dataX=get(ax_,'XData');
set(ax_,'XData',handles.Time,'YData',dataY); set(handles.axes1,'XLim',[handles.Time(1) handles.Time(end)]);
ax_=get(handles.axes2,'Children'); dataY=get(ax_,'YData'); dataX=get(ax_,'XData');
set(ax_,'XData',handles.Time,'YData',dataY); set(handles.axes1,'XLim',[handles.Time(1) handles.Time(end)]);
ax_=get(handles.axes3,'Children'); dataY=get(ax_,'YData'); dataX=get(ax_,'XData');
set(ax_,'XData',handles.Time,'YData',dataY); set(handles.axes1,'XLim',[handles.Time(1) handles.Time(end)]);

guidata(hObject, handles);
plotta_ora(handles)
end




%% --- Executes on button press in decimate.
function decimate_Callback(hObject, eventdata, handles)
% hObject    handle to decimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~any(strcmp(fieldnames(handles),'column')); err=msgbox('no data stored!'); waitfor(err); return; end

set(hObject,'Enable','on','BackGroundColor','white')
handles.Ndec=str2num(get(hObject,'String'))
%if any(strcmp(fieldnames(handles),'dec')); k=menu('decimate again?','yes','no');end
%if k==1
for n=1:length(handles.column)
    eval(['handles.' handles.column{n} ' = downsample(handles.' handles.column{n} ',' num2str(handles.Ndec) ');'])
end

%elseif k==2
%    return
%end

handles.dec='ok';

guidata(hObject, handles);
plotta_ora(handles);
end
%set(hObject,'String','decimate','BackgroundColor',[0.75 0.75 0.75]);

function decimate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to decimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[0.75 0.75 0.75]);
    waitfor(hObject,'String')
end
end

%% --- Executes on button press in smooth

function smooth_Callback(hObject, eventdata, handles)
% hObject    handle to smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of smooth as text
%        str2double(get(hObject,'String')) returns contents of smooth as a double

set(hObject,'Enable','on','BackGroundColor','white')
handles.sm=str2num(get(hObject,'String'));

if isempty(handles.sm); handles.sm=500; end
ginput(1)
ax_=get(gcf,'CurrentAxes');

for n=1:3
    eval(['s=find(ax_==handles.axes' num2str(n) ');'])
    if s==1; break; end
end
finestra=n;

posy=get(ax_,'Ylim');
posx=get(ax_,'Xlim');

eval(['h_=findobj(''Tag'',''edit' num2str(n) ''');']); %n=numero asse
s=str2double(get(h_,'String'));  %numero della colonna

%if any(strcmp(handles.column,{[handles.column{s} 'o']}));
%s=find(strcmp(handles.column,{[handles.column{s} 'o']}));
%end

%running mean
eval(['pippo=(handles.' handles.column{s} ');']);
if (handles.sm)/2==floor(handles.sm/2); handles.sm=handles.sm+1; end %check sul numero dispari
l=(handles.sm-1)/2; %es:(101-1)/2=50;

pm=pippo(l+1:length(pippo)-l);
for i=l+1:length(pippo)-l;
    pm(i-l)=sum(pippo(i-l:i+l));
end
pippo(l+1:length(pippo)-l)=pm/(handles.sm);

%windowSize = handles.sm;
%pm2=filter(ones(1,windowSize)/windowSize,1,pippo);
%smooths data Y using a handles.sm -point moving average.
%eval(['handles.' handles.column{s} '=smooth(handles.' handles.column{s} ',handles.sm);'])

eval(['handles.' handles.column{s} 'o=handles.' handles.column{s} ';'])
eval(['handles.' handles.column{s} '=pippo;'])
if ~strcmp(handles.column,{[handles.column{s} 'o']});
    handles.column(end+1)={[handles.column{s} 'o']};
end


%handles.g1=find(strcmp(handles.column,{[handles.column{s} 'smooth']}));

h_=findobj('Tag','edit1LB'); set(h_,'String',handles.column);
h_=findobj('Tag','edit2LB'); set(h_,'String',handles.column);
h_=findobj('Tag','edit3LB'); set(h_,'String',handles.column);



guidata(hObject, handles);

plotta_ora(handles);
end

% --- Executes during object creation, after setting all properties.
function smooth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

%% --- Executes on button press in cut_dt

function cut_dt_Callback(hObject, eventdata, handles)
% hObject    handle to cut_dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cut_dt as text
%        str2double(get(hObject,'String')) returns contents of cut_dt as a double
set(hObject,'Enable','on')


handles.dt=str2double(get(hObject,'String'));
%set(hObject,'String',handles.dt,'BackgroundColor',[0.75 0.75 0.75])
ll=find(handles.Stamp(:,1)==handles.dt); %,1,'first');
if length(ll) <= 100; h=msgbox('attention: number of residuals less than 100'); waitfor(h); return; end
%nn=find(handles.Stamp(:,1)==handles.dt,1,'last');

for n=1:length(handles.column)
    eval(['handles.' handles.column{n} ' = handles.' handles.column{n} '(ll,1);'])
end


% ax_=get(handles.axes1,'Children'); dataY=get(ax_,'YData'); dataX=get(ax_,'XData');
% set(ax_,'XData',dataX(ll),'YData',dataY(ll)); set(handles.axes1,'XLim',[dataX(ll(1)) dataX(ll(end))]);
% ax_=get(handles.axes2,'Children'); dataY=get(ax_,'YData'); dataX=get(ax_,'XData');
% set(ax_,'XData',dataX(ll),'YData',dataY(ll)); set(handles.axes2,'XLim',[dataX(ll(1)) dataX(ll(end))]);
% ax_=get(handles.axes3,'Children'); dataY=get(ax_,'YData'); dataX=get(ax_,'XData');
% set(ax_,'XData',dataX(ll),'YData',dataY(ll)); set(handles.axes3,'XLim',[dataX(ll(1)) dataX(ll(end))]);


set(hObject,'String','cut_dt')
guidata(hObject, handles);
plotta_ora(handles)
end

% --- Executes during object creation, after setting all properties.
function cut_dt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cut_dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



%% --- Executes on button press in cut.
function cut_Callback(hObject, eventdata, handles)
% hObject    handle to cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hOb=findobj('Tag','XLab');
h_ele=get(hOb,'Value');

if h_ele==2
    t_cut=handles.Time/1000;
elseif h_ele==3;
    if any(strcmp(fieldnames(handles),'slip'))
        t_cut=handles.slip;
    elseif any(strcmp(fieldnames(handles),'Slip_Enc_2'))
        t_cut=handles.Slip_Enc_2;
    end
elseif h_ele==1
    t_cut=handles.XLab;
end

%left to right clicking sequence on plot
[xi,yi]=ginput(2) ;


mat(1,:)=abs(t_cut-ones(size(t_cut))*xi(1));
mat(2,:)=abs(t_cut-ones(size(t_cut))*xi(2));
ll(:,1)=find(mat(1,:)==min(mat(1,:)));
ll(:,2)=find(mat(2,:)==min(mat(2,:)));

handles.cutted(1,1)=handles.RateZero(ll(1,1));
handles.cutted(1,2)=handles.RateZero(ll(1,2));

for n=1:length(handles.column)
    eval(['handles.' handles.column{n} ' =handles.' handles.column{n} '(ll(:,1):ll(:,2),:);' ])
end


guidata(hObject, handles);
plotta_ora(handles)
end



%% --- Executes on button press in fft.
function fft_Callback(hObject, eventdata, handles)
% hObject    handle to fft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hOb=findobj('Tag','XLab');
h_ele=get(hOb,'Value');

if h_ele==2
    t_cut=handles.Time/1000;
elseif h_ele==3;
    if any(strcmp(fieldnames(handles),'slip'))
        t_cut=handles.slip;
    elseif any(strcmp(fieldnames(handles),'Slip_Enc_2'))
        t_cut=handles.Slip_Enc_2;
    end
elseif h_ele==1
    t_cut=handles.XLab;
end


ax_=get(gcf,'CurrentAxes');

%h_=findobj('Tag','dt_value');
ginput(1)
for n=1:3
    eval(['s=find(ax_==handles.axes' num2str(n) ');'])
    if s==1; break; end
end
finestra=n;

ax_=get(gcf,'CurrentAxes');
posy=get(ax_,'Ylim');
posx=get(ax_,'Xlim');
I1=find(abs(t_cut-posx(1))==min(abs(t_cut-posx(1))));
I2=find(abs(t_cut-posx(2))==min(abs(t_cut-posx(2))));

eval(['h_=findobj(''Tag'',''edit' num2str(n) ''');']); %n=numero asse
s=str2double(get(h_,'String'));  %numero della colonna


eval(['mmed_torq=mean(smooth(handles.' handles.column{s} '(1:100,:)));']);
eval(['handles.' handles.column{s} '=handles.' handles.column{s} '- mmed_torq;']);

%tapering

eval(['handles.' handles.column{s} '=handles.' handles.column{s} '(I1:I2).*hamming(length(handles.Stamp(I1:I2)));'])
dt=mode(handles.Stamp(I1:I2)); %str2double(get(h_,'String'));

nfft=4096;
eval(['Y =fft(handles.' handles.column{s} ',' num2str(nfft) ');'])
f = 1000/dt*(0:nfft/2)/nfft;
Pyy = Y.* conj(Y) / nfft;

plot(f,Pyy(1:nfft/2+1),'Parent',handles.axes4);
set(handles.axes4,'XLim',[0 250])

plot(1./f*1000/dt,Pyy(1:nfft/2+1),'Parent',handles.axes5);
set(handles.axes5,'XLim',[0 500])

set(hObject,'Value',0)
end
%eval(['plot(handles.Timestamp,handles.' handles.column{s} ...
%    ',''ob'',''parent'',handles.axes' num2str(finestra) ');']);
%eval(['legend(handles.axes' num2str(finestra) ...
%    ',[handles.column{' num2str(s) '}])'])



%% SAVE--------------------------------------------------------------------
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pat0=pwd;
[nome,pat]=uiputfile( ...
    {'*.m;*.fig;*.mat;*.mdl', 'All MATLAB Files (*.m, *.fig, *.mat, *.mdl)'; ...
    '*.*',                   'All Files (*.*)'}, ...
    'Save as',[pat0 '/' handles.filename]);
cd (pat)

handles.save=handles.column;


for j=1:length(handles.save)
    if j==length(handles.save)
        M(j,1)={['''' handles.save{j} '''']};
    else
        M(j,1)={['''' handles.save{j} '''' ',']};
    end
end
M1=cell2mat(M');

%for j=1:length(handles.column)
%    if j==length(handles.column)
%        O(j,1)={['''v' num2str(j) '''']};
%    else
%        O(j,1)={['''v' num2str(j) '''' ',']};
%    end

%    O1=cell2mat(O');
%end

%file header
name4=['header', nome];
%fid1 = fopen(name4,'wt');

%fprintf(fid1,' load=%d\n shearT=%d\n trigger=%d\n cuttingfrom=%d to=%d\n dt=%f\n smooth=%d\n', ...
%    handles.loadT, handles.shearT,handles.triggered,handles.cutted(1,1),handles.cutted(1,2),handles.dt,handles.sm);
%fprintf(fid1,' load=%d\n shearT=%d\n trigger=%d\n cuttingfrom=%d to=%d\n dt=%f\n smooth=%d\n', ...
%    handles.loadT, handles.shearT,handles.triggered,handles.dt,handles.sm);
%fclose(fid1);

nome2=[nome, 'RED.mat'];
%nome3=['originali', nome];


eval(['save(nome2,''-struct'',''handles'',' M1 ');'])
end
%eval(['save(nome3,''-struct'',''handles'',' O1 ');'])
%save('parametri','-struct','handles','loadT','shearT','triggered','cutted'
%,'dt','sm')


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end




%% --------------------------------------------------------------------
function Interactive_Callback(hObject, eventdata, handles)
% hObject    handle to Interactive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
keyboard
return
end



%% --- Executes on button press in calibration.
function calibration_Callback(hObject, eventdata, handles)
% hObject    handle to calibration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%definisce la calibrazion

if any(strcmp(fieldnames(handles),'smooth'));
    %msgbox('I will work only on not smoothed data')
    nomi=handles.column;
    for i=1:length(nomi)
        if strfind(nomi{i},'smooth'); handles=rmfield(handles,nomi{i});
            K=find(~strcmp(nomi{i},handles.column));
            handles.column=handles.column(K);
        end
    end
end


if any(strcmp(fieldnames(handles),'new'))
    
    nomi=handles.new;
    for i=1:length(nomi)
        handles = rmfield(handles, nomi(i));
        K=find(~strcmp(nomi(i),handles.column));
        handles.column=handles.column(K);
    end
    if any(strcmp(fieldnames(handles),'new')); handles=rmfield(handles,'new');end
end



dint=findobj('Tag','Rint'); rint=str2double(get(dint,'String'))/2000;
dext=findobj('Tag','Rext'); rext=str2double(get(dext,'String'))/2000;

contents=get(hObject,'Value');

%NB in contents 2-4 cal.e2 is an average between the two corrections for internal and external
%velocity input

%remember that handles.slope and handles.offset now contain the calibrations
%read from the .csv file and are referred to the measurements in the order:
% {'Time','VerticalLoad','VerticalDisplacement','ShearTorque','Angle','Tachometer','PulsCounter','DigitalFunctionGenerator'}

if contents==2; cal.t=1; cal.a=1; cal.l=1; cal.e=1.49997; cal.e2=0.998; cal.ag=1; cal.fg=1;  end %corrections
if contents==3; cal.t=1; cal.a=1; cal.l=1; cal.e=0.99997; cal.e2=0.998; cal.ag=1; cal.fg=1;  end %corrections
if contents==4; cal.t=1; cal.a=1; cal.l=1; cal.e=1.030897; cal.e2=1; cal.ag=1; cal.fg=1; end %corrections
if contents==5; cal.t=1; cal.a=1; cal.l=1; cal.e=1; cal.e2=1; cal.ag=1; cal.fg=1; end %calibration in DCS
if contents==6; cal.t=0.1287; cal.a=2.5e-3; cal.l=1e-3; cal.e=99.997; cal.e2=199.83; cal.ag=35.5; cal.fg=1; end % calibration removed

% calculations

cal.torque(1:10)=cal.t*3/2/pi/(rext^3-rint^3)*1E-6;
cal.ax(1:10)=cal.a*(100*9.807)/(pi*(rext^2-rint^2))*1E-6;
cal.lv(1:5)=cal.l*-0.001; cal.lv(6:10)=-cal.l;
cal.enc(1:10)=cal.e*(4/3*pi*((rext^2+rint*rext+rint^2)/(rext+rint))); 
cal.enc2(1:10)=cal.e2*(4/3*pi*((rext^2+rint*rext+rint^2)/(rext+rint))/60);
cal.ang(1:10)=cal.ag*(4/3*pi*((rext^2+rint*rext+rint^2)/(rext+rint))); 

%% ELABORAZIONE DATI

% Normal stress
b=(strfind(handles.column,'Load')); j=0; n=[];
for i=1:length(b); if ~isempty(b{i}); j=j+1; n(j)=i; end; end
for j=1:length(n)
    new.Normal=handles.(handles.column{n(j)})*cal.ax(j);
end

% Shear, Mu
b=(strfind(handles.column,'Torque')); j=0; n=[];
for i=1:length(b); if ~isempty(b{i}); j=j+1; n(j)=i; end; end
for j=1:length(n)
    new.Shear=handles.(handles.column{n(j)})*cal.torque(j);
    new.Mu=new.Shear./new.Normal;
end

% Vertical displacement-thickness (if zeroed)
b=(strfind(handles.column,'VerticalDisplacement')); j=0; n=[];
for i=1:length(b); if ~isempty(b{i}); j=j+1; n(j)=i; end; end
for j=1:length(n)
    new.Vdisp=handles.(handles.column{n(j)})*cal.lv(j);
end

% Thickness

h_=findobj('Tag','thick');
statoTHICK=get(h_,'Value');

h_=findobj('Tag','thickz');
thickz=get(h_,'String');
thickz=str2double(thickz);

if statoTHICK==1
    
    b=(strfind(handles.column,'VerticalDisplacement')); j=0; n=[];
    for i=1:length(b); if ~isempty(b{i}); j=j+1; n(j)=i; end; end
    for j=1
        new.Thickness=((handles.(handles.column{n(j)})-thickz)*cal.l);
    end
    
else %keyboard
end

%% Slip: Encoder OR Angle depends on selection RIFARE

% b=(strfind(handles.column,'Angle')); j=0; n=[];
% for i=1:length(b)
%     if ~isempty(b{i}); j=j+1; n(j)=i;
%     end
% end
% for j=1:length(n)
%     new.SlipA=handles.(handles.column{n(j)})*cal.ang(j)/360; %slip in m
% end

b=(strfind(handles.column,'Counter')); j=0; n=[];
for i=1:length(b)
    if ~isempty(b{i}); j=j+1; n(j)=i;
    end
end
for j=1:length(n)
    new.SlipPC=handles.(handles.column{n(j)})*cal.enc(j); %slip in m
end

% True angle ????
b=(strfind(handles.column,'Angle')); j=0; n=[];
for i=1:length(b); if ~isempty(b{i}); j=j+1; n(j)=i; end; end
for j=1:length(n)
    new.SlipAIncr=unwrap(handles.(handles.column{n(j)})*pi/180)*cal.ang(j); %slip in m, incremental
end

%Equivalent velocity from rpm
b=(strfind(handles.column,'Tachometer')); j=0; n=[];
for i=1:length(b); if ~isempty(b{i}); j=j+1; n(j)=i; end; end
for j=1:length(n)
    new.VelocityPC=handles.(handles.column{n(j)})*cal.enc2(j);
end

%% MISSING VELOCITY FROM SlipAIncr
addpath('C:\Users\Stefano\Dropbox\Stefano\Ricerca\RoSA\oct19calibrations\m2018')
time=handles.Time;
slip=new.SlipAIncr; 
slip=slip-slip(1,1);
slipsg=sgolayfilt(slip,3,101);
slipbf=brutal_filter_fx(31,1,100,slip,handles.Stamp);
% figure;plot(time,slip,'-r');hold on;plot(time,slipsg,'-g'); hold on; plot(time,slipbf,'-b');
% yyaxis right; plot(time(1:end-1),diff(slipsg)/handles.Stamp(1,1),'-g')
% yyaxis right; plot(time(1:end-1),diff(slipbf)/handles.Stamp(1,1),'-b')

new.VelocityAIncr=zeros(size(handles.Angle,1),size(handles.Angle,2));
new.VelocityAIncr(2:end,1)=diff(slipbf)/handles.Stamp(1,1);

% h=findobj('Tag','nodeEnc1');
% node=str2num(get(h,'String'));
% node1=node(:,1);
% f1crat=node(:,2);
%
% h=findobj('Tag','nodeEnc2');
% node=str2num(get(h,'String'));
% node2=node(:,1);
% f2crat=node(:,2);
%
% h=findobj('Tag','MaxEnc');
% maxE=str2double(get(h,'String')) ;
%
% b=(strfind(handles.column,'Encoder')); j=0; n=[];
% for i=1:length(b); if ~isempty(b{i}); j=j+1; n(j)=i; end; end
% for j=1:length(n);

%  clear d0 d1 d2 v vel
%  d0=handles.(handles.column{n(j)}) %d0 is the initial, raw encoder data
% Dd0 = diff(d0);
% Dd0(Dd0 == 0) = [];
% min(Dd0)
% if abs(min(Dd0)) < 1/4000
%     disp('Fine encoder')
% else
%     disp('Coarse encoder')
% end
% h_ele=handles.Stamp/1000;
% % qui provo a ridurre il problema del cambiamento in corsa di acq xlab.
% % per rirpistinare la versione di rosa2019 precedente imposta I_ele=[] e commenta la riga di seguito.
%
% hOb=findobj('Tag','AdjRate'); f_ele=get(hOb,'Value');
% I_ele=find(diff(h_ele)<mode(diff(h_ele)));
% testE=2;
%
% %=== se necessario aggiustare il rate
% if ~isempty(I_ele) && f_ele==1
%     testE=1; Fs=[];
%     sl_ele=d0;
%     I_range=[I_ele: I_ele+1];
%     I_fit=[I_ele-10:I_ele];
%     I_inter=[I_ele:I_ele+10];
%     cv=fit(handles.Time(I_fit)/1000, sl_ele(I_fit),'linear');
%     deltaslip=sl_ele(I_ele+1) - cv(handles.Time(I_ele+1)/1000);
%     sl_ele(I_ele+1:end)=sl_ele(I_ele+1:end)/1000 - deltaslip(1);
%     d0=sl_ele;
% end
%
%
% %che i massimi siano fuori dal rumore
% J=local_max(d0); I=[];
% if length(J) > 10
% [c,d]=hist(diff(d0(J))); Jc=find(c==max(c));
% r=abs(diff(d0(J))-d(Jc(1))); In=find(r > mode(diff(d))); I=J(In);
% end
%
% I=J; %cambiato 15/03
% %figure; plot(d0); hold on; plot(I,d0(I),'*r')
%
% fO=findobj('Tag','incremental'); fOV=get(fO,'value')
% if fOV==1
% if ~isempty(I)
% for ii=1:length(I)
%     if ii==length(I) & I(ii) + 2 < length(d0)
% d0(I(ii)+1:end)=d0(I(ii)+1:end)-d0(I(ii)+2)+d0(I(ii));
%     elseif ii==length(I) & I(ii) +2 >= length(d0)
%         break
%     else
% d0(I(ii)+1:I(ii+1))=d0(I(ii)+1:I(ii+1))-d0(I(ii)+2)+d0(I(ii));
%     end
% end
% end
% end %if fOV
%
% d0(d0<0)=0;
%
%     ITorq=findobj('Tag','Torque');
%     ITorqO=get(ITorq,'Value');
%
% h=handles.Stamp/1000;
% d1=d0*cal.enc(1);
% figure(22)
% plot(d1)
% if ITorqO==0
%     fcamp=(1./max(1/fref,h));
%     Fs=max(fcamp);
%     alphas=ceil(log10(Fs));
%     fc=max(fref/100,Fs/10^(alphas));
%     disp(['fc= ',num2str(fc)])
%     nodes=node2; fcrat=f2crat;
%     if i==1; nodes=node1; fcrat=f1crat;
%     end
%     if Fs > fref/100
%     disp('red')
%     %riduci l'effetto window
%     d=fdesign.lowpass('N,Fc',nodes,fcrat,Fs);
%
%     Hd=design(d,'window','Window',tukeywin(nodes+1,0));
%     d_sm=filtfilt(Hd.Numerator,1,d1);
%     else
%     d_sm=d1;
%     end
% else
% d_sm=d1;
% end
%
%
%
% % na=1;
% % dtn=mode(h)./na;
% %
% % if na > 1
% %     % serve a ricampionare e calcolare la V con maggiore dettaglio.
% %     % Cambiare na all'occorrenza
% % xx=handles.Time(1)/1000:dtn:handles.Time(end)/1000;
% % cs=spline(handles.Time/1000, d_sm);
% % d_sm1=ppval(cs,xx);
% % bb=zeros(size(d_sm1));
% % bb(2:end)=diff(d_sm1); bb(1)=bb(2);
% % h_1=diff(xx); h1=h_1; h1(2:end+1)=h_1;
% % v1=bb./h1;
% % v=resample(v1,1,na)';
% % d_sm=resample(d_sm1, 1, na)';
% % a=length(v) - length(h);
% % for ia=1:a
% % v(end+1)=v(end);
% % d_sm(end+1)=d_sm(end);
% % end
% % else
% bb=zeros(size(d_sm));
% bb=diff(d_sm); bb(end+1)=bb(end);
% v=(bb./h)./handles.tconv;
% % end
%
% eval(['new.SlipVel_Enc_' num2str(j) '=v;']);
% eval(['new.Slip_Enc_' num2str(j) '=(d_sm);']);
%
%
% %if j==2; v=smooth(v,300); end
% if j==1; del0=d_sm; v0=v; end
% if j==2;
% del10=max(d_sm,del0);
% Iv=find(v < 0.01 ); %% cambiare
% d_sm(Iv)=del0(Iv);
% v=max(v0,v);
% del=max(del0,d_sm);
% %if j==2; adel=abs(d_sm-del0)./del0.*100; Iadel=find(adel < 1);
% %del(Iadel)=del0(Iadel);
% %end
% %bb=zeros(size(del));
% %bb(2:end)=diff(del); bb(1)=bb(2);
%
% %if j==1; Ivmax=find(v<=0.05); new.vel=zeros(size(v)); new.vel(Ivmax)=v(Ivmax); end
% %if j==2; Ivmax=find(v>0.05);  new.vel(Ivmax)=new.SlipVel_Enc_2(Ivmax); end
% new.slip=del10;
% if max(v) <= 20E-3
%     new.vel=v0;
% else
% new.vel=v; %bb./h;
% end
%
% %find outlayers in vel
% IOL=find(new.vel >=10);
% for i=1:length(IOL)
% new.vel(IOL(i))=new.vel(IOL(i)-1);
% end
% IOL=find(new.vel < -1);
% for i=1:length(IOL)
% new.vel(IOL(i))=new.vel(IOL(i)-1);
% end
%
% end
% end
%
%% temp estimation
% dn=50;
% time2=cumsum(handles.Stamp);
%
% if any(strcmp(fieldnames(new),'vel'))
%     [Temp]=temp(handles.Time/1000,new.vel,new.shear1,new.slip, dn); %change for vel
%
% elseif max(new.vel) <= 20e-3
%     [Temp]=temp(handles.Time/1000,new.SlipVel_Enc_1,new.shear1,new.Slip_Enc_1,dn);
% else
%     [Temp]=temp(handles.Time/1000,new.SlipVel_Enc_2,new.shear1,new.Slip_Enc_2,dn); %change for vel
% end
% new.TempE=interp1(time2(1:dn:end),Temp,time2);

%% write in handles
nomi=[];

[a,b]=size(handles.column); [aa,bb]=size(fieldnames(new));
inp1=handles.column;
if aa==1 && aa==b || bb==1 && a==1; inp1=handles.column'; end
nomi=[inp1 ; fieldnames(new)];

handles.column=[];
handles.column=nomi';

nomi2=fieldnames(new);
for i=1:length(nomi2)
    eval(['handles.' char(nomi2(i)) '=new.' char(nomi2(i)) ';'])
end

h_=findobj('Tag','edit1LB'); set(h_,'String',handles.column)
h_=findobj('Tag','edit2LB'); set(h_,'String',handles.column)
h_=findobj('Tag','edit3LB'); set(h_,'String',handles.column)

handles.new=fieldnames(new)';

guidata(hObject, handles);

end


%% definisce il diametro interno ed esterno
function Rint_Callback(hObject, eventdata, handles)
% hObject    handle to Rint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% Hints: get(hObject,'String') returns contents of Rint as text
%        str2double(get(hObject,'String')) returns contents of Rint as a double

% --- Executes during object creation, after setting all properties.
function Rint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function Rext_Callback(hObject, eventdata, handles)
% hObject    handle to Rext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Rext as text
%        str2double(get(hObject,'String')) returns contents of Rext as a double
end


% --- Executes during object creation, after setting all properties.
function Rext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



%% --- Executes on button press in Gefran.
function Gefran_Callback(hObject, eventdata, handles)
% hObject    handle to Gefran (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Gefran
stato=get(hObject,'Value');

if stato==1
    ButtonName = questdlg('GEFRAN preliminary operations: (1) identify the main sampling rate with cut_dt (2) decimate the data to 1s', ...
        'Do you want to proceed now?', 'No', 'Yes','No');
    
    switch ButtonName,
        case 'No',
            disp('no');
            set(hObject,'Value',0);
        case 'Yes',
            
            I=find(~strcmp(handles.column,'SpeedGEF')); handles.column=handles.column(I);
            I=find(~strcmp(handles.column,'timeGEF')); handles.column=handles.column(I);
            I=find(~strcmp(handles.column,'TorqueGEF')); handles.column=handles.column(I);
            
            Path2='/media/disk1/shivadir/Shiva Experiments';
            name=handles.filename(1:4);
            
            list=dir([Path2 '/' name '*']);
            if ~isempty(list); Path2=[Path2 '/' list.name '/'];
                name2=dir([Path2 '*.txt']);
                FileGEF=name2.name;
            else
                [FileGEF,Path2] = uigetfile( '/media/disk1/shivadir/*.txt', ...
                    'Multiple File Detected: Select the file GEFRAN to load');
                name2='';
            end
            
            
            gefran1=1; k=0 ;
            while ~isstruct(gefran1)
                k=k+1;
                gefran1=importdata([Path2 '/' FileGEF],'\t',k);
            end
            
            I=find(strncmp(gefran1.textdata,'Time	Speed',6)); if ~isempty(I); new.timeGEF=gefran1.data(:,1); new.VGEF(:,1)=gefran1.data(:,2); end
            I=find(strncmp(gefran1.textdata,'Act Torque',6)); if ~isempty(I); new.TqGEF(:,1)=gefran1.data(:,I(1)-1); end
            I=find(strncmp(gefran1.textdata,'Speed',5)); if ~isempty(I); new.VGEF(:,1)=gefran1.data(:,I(1)-1); end
            I=find(strncmp(gefran1.textdata,'time',4)); if ~isempty(I); new.timeGEF(:,1)=gefran1.data(:,I(1)-1); end
            
            dt=diff(new.timeGEF); dt(end+1)=dt(1);
            
            %for j=1:length(handles.column);
            %    eval(['new.' handles.column{j} '= downsample(handles.' handles.column{j} ',' num2str(25) ');'])
            %end
            %dati di calibrazione
            
            nomi=[];
            [a,b]=size(handles.column); [aa,bb]=size(fieldnames(new));
            inp1=handles.column;
            if aa==1 & aa==b | bb==1 & a==1; inp1=handles.column'; end
            nomi=[inp1 ; fieldnames(new)];
            
            handles.column=[];
            handles.column=nomi';
            nomi2=fieldnames(new);
            
            for k=1:length(nomi2)
                eval(['handles.' char(nomi2(k)) '=new.' char(nomi2(k)) ';'])
            end
            
            guidata(hObject, handles);
            
    end
end % switch
end

%%%%%%%%%%% fine GEF


%% --- Executes on button press in running.
function running_Callback(hObject, eventdata, handles)
% hObject    handle to running (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hOb=findobj('Tag','XLab');
h_ele=get(hOb,'Value');

if h_ele==2
    t_cut=handles.Time/1000;
elseif h_ele==3;
    if any(strcmp(fieldnames(handles),'slip'))
        t_cut=handles.slip;
    elseif any(strcmp(fieldnames(handles),'Slip_Enc_2'))
        t_cut=handles.Slip_Enc_2;
    end
elseif h_ele==1
    t_cut=handles.XLab;
end



[xi,yi]=ginput(2) ;
mat(1,:)=abs(t_cut-ones(size(t_cut))*xi(1));
mat(2,:)=abs(t_cut-ones(size(t_cut))*xi(2));
ll(:,1)=find(mat(1,:)==min(mat(1,:)));
ll(:,2)=find(mat(2,:)==min(mat(2,:)));


ax_=get(gcf,'CurrentAxes');
for n=1:3
    eval(['s=find(ax_==handles.axes' num2str(n) ');'])
    if s==1; break; end
end
finestra=n;
posy=get(ax_,'Ylim');
posx=get(ax_,'Xlim');

eval(['h_=findobj(''Tag'',''edit' num2str(n) ''');']); %n=numero asse
s=str2double(get(h_,'String'));                        %numero della colonna

%cla

eval(['mmed_torq=mean(smooth(handles.' handles.column{s} '(ll(:,1):ll(:,2),:)));']);
eval(['handles.' handles.column{s} '=handles. ' handles.column{s} '- mmed_torq;']);
eval(['plot(handles.XLab,handles.' handles.column{s} ',''ob'',''parent'',handles.axes' num2str(finestra) ''');']);
set(ax_,'Xlim',posx)

guidata(hObject, handles);
plotta_ora(handles)
end


% --- Executes on button press in fluid.
function fluid_Callback(hObject, eventdata, handles)
% hObject    handle to fluid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fluid
end


% --------------------------------------------------------------------
function print_Callback(hObject, eventdata, handles)
% hObject    handle to file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hf=figure;
hn_=copyobj(handles.axes1,hf); h1=get(hn_,'Position'); h1(1)=h1(1)+0.1; set(hn_,'Position',h1);
ah_=get(handles.axes1,'children'); nom=get(ah_,'DisplayName');
hnl_=get(hn_,'YLabel'); set(hnl_,'string',nom)

hn_=copyobj(handles.axes2,hf); h1=get(hn_,'Position'); h1(1)=h1(1)+0.1; set(hn_,'Position',h1);
ah_=get(handles.axes2,'children'); nom=get(ah_,'DisplayName');
hnl_=get(hn_,'YLabel'); set(hnl_,'string',nom)

hn_=copyobj(handles.axes3,hf); h1=get(hn_,'Position'); h1(1)=h1(1)+0.1; set(hn_,'Position',h1);
ah_=get(handles.axes3,'children'); nom=get(ah_,'DisplayName');
hnl_=get(hn_,'YLabel'); set(hnl_,'string',nom)
end

%ax_=get(handles.axes1,'Parent');
%ah_=get(ax_,'Children')
%hf=figure;
%for n=1:length(ah_)
%eval(['h_new = copyobj(ah_(' num2str(n) '),hf);']);
%end




% --- Executes on button press in outlayers.
function outlayers_Callback(hObject, eventdata, handles)
% hObject    handle to outlayers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%zoom xon;
hOb=findobj('Tag','XLab');
h_ele=get(hOb,'Value');

if h_ele==2
    t_cut=handles.Time/1000;
elseif h_ele==3;
    if any(strcmp(fieldnames(handles),'slip'))
        t_cut=handles.slip;
    elseif any(strcmp(fieldnames(handles),'Slip_Enc_2'))
        t_cut=handles.Slip_Enc_2;
    end
elseif h_ele==1
    t_cut=handles.XLab;
end

button=1; i=0;
while button==1
    i=i+1
    [xi(i),yi,button]=ginput(1) ;
    button
end


ax_=get(gcf,'CurrentAxes');
for n=1:3
    eval(['s=find(ax_==handles.axes' num2str(n) ');'])
    if s==1; break; end
end
finestra=n;
posy=get(ax_,'Ylim');
posx=get(ax_,'Xlim');

eval(['h_=findobj(''Tag'',''edit' num2str(n) ''');']); %n=numero asse
s=str2double(get(h_,'String'));  %numero della colonna


for i=1:length(xi);
    mat(1,:)=abs(t_cut-ones(size(t_cut))*xi(i));
    ll(i)=find(mat(1,:)==min(mat(1,:)),1,'first');
end

if button==3
    
    for i=1:length(xi);
        %running mean
        
        if ll==1; eval(['handles.' handles.column{s} '(ll(i))=handles.' handles.column{s} '(ll(i)+1);']);
        else
            eval(['handles.' handles.column{s} '(ll(i))=handles.' handles.column{s} '(ll(i)-1);']);
        end
        
    end
    
    
elseif button==2
    eval(['handles.' handles.column{s} '(ll(1):ll(end))=handles.' handles.column{s} '(ll(1)-1);']);
    
end %if button


guidata(hObject, handles);

plotta_ora(handles)
end




% --- Executes on selection change in edit1LB.
function edit1LB_Callback(hObject, eventdata, handles)
% hObject    handle to edit1LB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns edit1LB contents as cell array
%        contents{get(hObject,'Value')} returns selected item from edit1LB


s=get(hObject,'Value');
h=findobj('Tag','edit1');
handles.g1=s;
set(h,'String',s);

guidata(hObject, handles);

plotta_ora(handles)
end



% --- Executes during object creation, after setting all properties.
function edit1LB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1LB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end




% --- Executes on selection change in edit2LB.
function edit2LB_Callback(hObject, eventdata, handles)
% hObject    handle to edit2LB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns edit2LB contents as cell array
%        contents{get(hObject,'Value')} returns selected item from edit2LB

s=get(hObject,'Value');
h=findobj('Tag','edit2');
handles.g2=s;
set(h,'String',s);

guidata(hObject, handles);

plotta_ora(handles)
end


% --- Executes during object creation, after setting all properties.
function edit2LB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2LB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on selection change in edit3LB.
function edit3LB_Callback(hObject, eventdata, handles)
% hObject    handle to edit3LB (see GCBO)
%
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns edit3LB contents as cell array
%        contents{get(hObject,'Value')} returns selected item from edit3LB

s=get(hObject,'Value');
h=findobj('Tag','edit3');
handles.g3=s;
set(h,'String',s);

guidata(hObject, handles);

plotta_ora(handles)
end



% --- Executes during object creation, after setting all properties.
function edit3LB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3LB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end






% --------------------------------------------------------------------
function Figure_Callback(hObject, eventdata, handles)
% hObject    handle to Figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

%% LOAD --------------------------------------------------------------------
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


ah_=get(handles.axes1,'children'); set(ah_,'XData',[],'YData',[]);
ah_=get(handles.axes2,'children'); set(ah_,'XData',[],'YData',[]);
ah_=get(handles.axes3,'children'); set(ah_,'XData',[],'YData',[]);

I=strcmp(fieldnames(handles),'column');

if any(I)
    for i=1:length(handles.column)
        eval(['handles=rmfield(handles,''', handles.column{i}, ''');'])
    end
    handles=rmfield(handles,'column');
    
end
clear file1
I=strcmp(fieldnames(handles),'new'); if any(I); handles=rmfield(handles,'new'); end
I=strcmp(fieldnames(handles),'X'); if any(I); handles=rmfield(handles,'X'); end
I=strcmp(fieldnames(handles),'TimeZero'); if any(I); handles=rmfield(handles,'TimeZero'); end



%definisce i grafici da plottare:
%qui ci sono i default
handles.g1=2;
handles.g2=3;
handles.g3=5;

ax_=findobj('Tag','edit1'); set(ax_,'String',handles.g1);
ax_=findobj('Tag','edit2'); set(ax_,'String',handles.g2);
ax_=findobj('Tag','edit3'); set(ax_,'String',handles.g3);

h_=findobj('Tag','edit1LB'); set(h_,'String',1);
h_=findobj('Tag','edit2LB'); set(h_,'String',1);
h_=findobj('Tag','edit3LB'); set(h_,'String',1);




[FileName,PathName] = uigetfile('*.*','All Files (*.*)', ...
    'C:\Users\Stefano\Dropbox\Ricerca\SHIVA');

cd (PathName)

data=load(FileName);


dataName=fieldnames(data);
% if length(dataName); data=getfield(data,dataName{1}); end


h_=findobj('Tag','dt_value');
stato=get(h_,'Value');
if isempty(stato)
    handles.dt=0.04;
else
    handles.dt=stato;
end


%set(h_,'String',handles.dt);

handles.filename=FileName;

handles.sm=0;
handles.triggered=0;
handles.cutted=[0 0];
handles.loadT=0;
handles.shearT=0;
ll=1;
nn=length(data.Time);

handles.column=fieldnames(data)

for i=1:length(handles.column)
    handles.(handles.column{i})=data.(handles.column{i}); %debuggato
    %     eval(['handles.' handles.column{i} '=data.' handles.column{i} ';'])
end


h_=findobj('Tag','edit1LB'); set(h_,'String',handles.column);
h_=findobj('Tag','edit2LB'); set(h_,'String',handles.column);
h_=findobj('Tag','edit3LB'); set(h_,'String',handles.column);

handles.load=1;
handles.Done=1;
new=handles;
guidata(hObject, handles);
guidata(hObject, new);

handles.zoom=0;

plotta_ora(handles);
end






%% BINARY
% --------------------------------------------------------------------
function binary_Callback(hObject, eventdata, handles)
% hObject    handle to binary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[nome,pat]=uiputfile( ...
    {'*.txt', 'All MATLAB Files (*txt)'; ...
    '*.*',                   'All Files (*.*)'}, ...
    'Write as',['~/',handles.filename]);
cd (pat)


for j=1:length(handles.column)
    C(j,1)={'%10.6f '};
    if j==length(handles.column)
        C(j,1)={'%10.6f\n'};
    end
end
C1=cell2mat(C');

for j=1:length(handles.column)
    N(j,1)={['handles.' handles.column{j} '(l,1),']};
    if j==length(handles.column)
        N(j,1)={['handles.' handles.column{j} '(l,1)']};
    end
end
N1=cell2mat(N');


for j=1:length(handles.column)
    M(j,1)={['''' handles.column{j} '''' ',']};
    if j==length(handles.column)
        M(j,1)={['''' handles.column{j} '''']};
    end
end
M1=cell2mat(M');

%for j=1:length(handles.column)
%    O(j,1)={['''v' num2str(j) '''' ',']};
%    if j==length(handles.column)
%         O(j,1)={['''v' num2str(j) '''']};
%    end
%end
%O1=cell2mat(O');


for j=1:length(handles.column)
    S(j,1)={'%s '};
    if j==length(handles.column)
        S(j,1)={'%s\n'};
    end
end
S1=cell2mat(S');

%write in a file
nome2=[nome, 'RED.txt'];
fid = fopen(nome2,'wt');
eval(['fprintf(fid,''' S1 ''',' M1 ');'])

eval(['len=length(handles.' handles.column{1} ');'])

for l=1:len
    eval(['fprintf(fid,''' C1 ''',' N1 ');'])
end
fclose(fid);
if ~ strcmp(fieldnames(handles),'dt'); msgbox(['ATTENTION: handles.dt=none']); end
end



% --- Executes on button press in slipON.
function slipON_Callback(hObject, eventdata, handles)
% hObject    handle to slipON (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of slipON
end



% --- Executes on button press in GH.
function GH_Callback(hObject, eventdata, handles)
% hObject    handle to GH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GH

end


% --- Executes on button press in TC.
function TC_Callback(hObject, eventdata, handles)
% hObject    handle to TC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TC
end



%% SAVERED--------------------------------------------------------------------
function saveRED_Callback(hObject, eventdata, handles)
% hObject    handle to saveRED (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pat0=pwd;
[nome,pat]=uiputfile( ...
    {'*.m;*.fig;*.mat;*.mdl', 'All MATLAB Files (*.m, *.fig, *.mat, *.mdl)'; ...
    '*.*',                   'All Files (*.*)'}, ...
    'Save as',[pat0 '/' handles.filename]);
cd (pat)


h_=findobj('Tag','fluid');
statoF=get(h_,'Value');
h_=findobj('Tag','GH');
statoGH=get(h_,'Value');

if statoF==1
    handles.save={'Time' 'shear1' 'EffPressure' 'Mu1' 'Pf' 'LVDT_low' 'LVDT_high' 'vel' 'slip' 'TempE' 'TempM'};
elseif statoGH==1
    handles.save={'Time' 'shear1' 'Normal' 'Mu1' 'dspring' 'LVDT_low' 'vel' 'slip'}; %'TempE' 'TempM'};
else
    handles.save={'Time' 'shear1' 'Normal' 'Mu1' 'LVDT_low' 'vel' 'slip' 'TempE'};
end


for j=1:length(handles.save)
    if j==length(handles.save)
        M(j,1)={['''' handles.save{j} '''']};
    else
        M(j,1)={['''' handles.save{j} '''' ',']};
    end
end
M1=cell2mat(M');

%for j=1:length(handles.column)
%    if j==length(handles.column)
%        O(j,1)={['''v' num2str(j) '''']};
%    else
%        O(j,1)={['''v' num2str(j) '''' ',']};
%    end

%    O1=cell2mat(O');
%end

%file header
name4=['header', nome];
%fid1 = fopen(name4,'wt');

%fprintf(fid1,' load=%d\n shearT=%d\n trigger=%d\n cuttingfrom=%d to=%d\n dt=%f\n smooth=%d\n', ...
%    handles.loadT, handles.shearT,handles.triggered,handles.cutted(1,1),handles.cutted(1,2),handles.dt,handles.sm);
%fprintf(fid1,' load=%d\n shearT=%d\n trigger=%d\n cuttingfrom=%d to=%d\n dt=%f\n smooth=%d\n', ...
%    handles.loadT, handles.shearT,handles.triggered,handles.dt,handles.sm);
%fclose(fid1);

nome2=[nome, '.mat'];
%nome3=['originali', nome];


eval(['save(nome2,''-struct'',''handles'',' M1 ');'])
%eval(['save(nome3,''-struct'',''handles'',' O1 ');'])
%save('parametri','-struct','handles','loadT','shearT','triggered','cutted'
%,'dt','sm')
end



function nodeEnc1_Callback(hObject, eventdata, handles)
% hObject    handle to nodeEnc1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nodeEnc1 as text
%        str2double(get(hObject,'String')) returns contents of nodeEnc1 as a double
end

% --- Executes during object creation, after setting all properties.
function nodeEnc1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nodeEnc1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function MaxEnc_Callback(hObject, eventdata, handles)
% hObject    handle to MaxEnc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxEnc as text
%        str2double(get(hObject,'String')) returns contents of MaxEnc as a double
end

% --- Executes during object creation, after setting all properties.
function MaxEnc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxEnc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function nodeEnc2_Callback(hObject, eventdata, handles)
% hObject    handle to nodeEnc2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nodeEnc2 as text
%        str2double(get(hObject,'String')) returns contents of nodeEnc2 as a double
end

% --- Executes during object creation, after setting all properties.
function nodeEnc2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nodeEnc2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function T0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function T0_Callback(hObject, eventdata, handles)
% hObject    handle to TT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TT as text
%        str2double(get(hObject,'String')) returns contents of TT as a double
end

% --- Executes during object creation, after setting all properties.
function TT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in AdjRate.
function AdjRate_Callback(hObject, eventdata, handles)
% hObject    handle to AdjRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AdjRate
end

% --- Executes on button press in Torque.
function Torque_Callback(hObject, eventdata, handles)
% hObject    handle to Torque (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Torque
end

% --- Executes on button press in off_enc_0.
function off_enc_0_Callback(hObject, eventdata, handles)
% hObject    handle to off_enc_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


I=handles.triggered;
handles.Encoder2(1:I)=0;
handles.Encoder(1:I)=0;
guidata(hObject, handles);
end

% --- Executes on button press in incremental.
function incremental_Callback(hObject, eventdata, handles)
% hObject    handle to incremental (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --- Executes on button press in vac.
function vac_Callback(hObject, eventdata, handles)
% hObject    handle to vac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vac
end

% --- Executes on button press in thick.
function thick_Callback(hObject, eventdata, handles)
% hObject    handle to thick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of thick
end


function thickz_Callback(hObject, eventdata, handles)
% hObject    handle to thickz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thickz as text
%        str2double(get(hObject,'String')) returns contents of thickz as a double
end

% --- Executes during object creation, after setting all properties.
function thickz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thickz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



% hOb=findobj('Tag','XLab');
% h_ele=get(hOb,'Value');
%
% if h_ele==2
%    handles.X=handles.Time/1000;
% elseif h_ele==3;
%         if any(strcmp(fieldnames(handles),'slip'))
%          handles.X=handles.slip;
%         elseif any(strcmp(fieldnames(handles),'Slip_Enc_2'))
%          handles.X=handles.Slip_Enc_2;
%         end
% elseif h_ele==1
%          handles.X=handles.Rate;
% end
% eval(['plot(handles.X,handles.' handles.column{(handles.g3)} ',''ob'',''parent'',handles.axes3);']);
% legend(handles.axes3,[handles.column{handles.g3}])

% --- Executes during object creation, after setting all properties.


% --- Executes on selection change in popupAI6.
function popupAI6_Callback(hObject, eventdata, handles)
% hObject    handle to popupAI6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupAI6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupAI6
end

% --- Executes during object creation, after setting all properties.
function popupAI6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupAI6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in popupAI7.
function popupAI7_Callback(hObject, eventdata, handles)
% hObject    handle to popupAI7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupAI7 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupAI7
end

% --- Executes during object creation, after setting all properties.
function popupAI7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupAI7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in popupAI8.
function popupAI8_Callback(hObject, eventdata, handles)
% hObject    handle to popupAI8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupAI8 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupAI8
end

% --- Executes during object creation, after setting all properties.
function popupAI8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupAI8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in popupAI9.
function popupAI9_Callback(hObject, eventdata, handles)
% hObject    handle to popupAI9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupAI9 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupAI9
end

% --- Executes during object creation, after setting all properties.
function popupAI9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupAI9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in popupAI10.
function popupAI10_Callback(hObject, eventdata, handles)
% hObject    handle to popupAI10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupAI10 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupAI10
end

% --- Executes during object creation, after setting all properties.
function popupAI10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupAI10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in popupAI16.
function popupAI16_Callback(hObject, eventdata, handles)
% hObject    handle to popupAI16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupAI16 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupAI16
end

% --- Executes during object creation, after setting all properties.
function popupAI16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupAI16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in popupAI17.
function popupAI17_Callback(hObject, eventdata, handles)
% hObject    handle to popupAI17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupAI17 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupAI17
end

% --- Executes during object creation, after setting all properties.
function popupAI17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupAI17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in popupPF.
function popupPF_Callback(hObject, eventdata, handles)
% hObject    handle to popupPF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupPF contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupPF
end

% --- Executes during object creation, after setting all properties.
function popupPF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupPF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in popupAI18.
function popupAI18_Callback(hObject, eventdata, handles)
% hObject    handle to popupAI18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupAI18 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupAI18
end

% --- Executes during object creation, after setting all properties.
function popupAI18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupAI18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in popupPC.
function popupPC_Callback(hObject, eventdata, handles)
% hObject    handle to popupPC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupPC contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupPC
end

% --- Executes during object creation, after setting all properties.
function popupPC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupPC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes when uipanel1 is resized.
function uipanel1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end
