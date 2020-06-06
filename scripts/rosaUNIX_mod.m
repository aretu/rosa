%% initialization code
function varargout = rosaUNIX_mod(varargin)
% ROSAUNIX_MOD M-file for rosaUNIX_mod.fig
%      ROSAUNIX_MOD, by itself, creates a new ROSAUNIX_MOD or raises the existing
%      singleton*.
%
%      H = ROSAUNIX_MOD returns the handle to a new ROSAUNIX_MOD or the handle to
%      the existing singleton*.
%
%      ROSAUNIX_MOD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ROSAUNIX_MOD.M with the given input arguments.
%
%      ROSAUNIX_MOD('Property','Value',...) creates a new ROSAUNIX_MOD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before shivaWIN_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to rosaUNIX_mod_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help rosaUNIX_mod

% Last Modified by GUIDE v2.5 21-Jan-2015 11:30:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @rosaUNIX_mod_OpeningFcn, ...
                   'gui_OutputFcn',  @rosaUNIX_mod_OutputFcn, ...
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

%% before GUI is visible

function rosaUNIX_mod_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to rosaUNIX_mod (see VARARGIN)

% Choose default command line output for rosaUNIX_mod
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using rosaUNIX_mod.

% UIWAIT makes rosaUNIX_mod wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = rosaUNIX_mod_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%% pushbutton1 = ?
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(~, ~, handles)
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


%% WRITE (scrive come txt i dati elaborati)
function write_Callback(~, ~, handles)
% hObject    handle to write (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[nome,pat]=uiputfile( ...
      {'*.txt', 'All MATLAB Files (*txt)'; ...
        '*.*',                   'All Files (*.*)'}, ...
         'Write as',['~/',handles.filename]);
cd (pat)


Iq=strcmp(handles.column,'vel'); if ~any(Iq); handles.vel=zeros(size(handles.Time)); 
handles.column{end+1}='vel'; end
Iq=strcmp(handles.column,'slip'); if ~any(Iq); handles.slip=zeros(size(handles.Time)); 
handles.column{end+1}='slip'; end
Iq=strcmp(handles.column,'SlipVel_Enc_2'); if ~any(Iq); handles.SlipVel_Enc_2=zeros(size(handles.Time)); 
handles.column{end+1}='SlipVel_Enc_2'; end
Iq=strcmp(handles.column,'Slip_Enc_2'); if ~any(Iq); handles.Slip_Enc_2=zeros(size(handles.Time)); 
handles.column{end+1}='Slip_Enc_2'; end
Iq=strcmp(handles.column,'LVDT_1'); if ~any(Iq); handles.LVDT_1=zeros(size(handles.Time)); 
handles.column{end+1}='LVDT_1'; end
Iq=strcmp(handles.column,'LVDT_2'); if ~any(Iq); handles.LVDT_2=zeros(size(handles.Time)); 
handles.column{end+1}='LVDT_2'; end
Iq=strcmp(handles.column,'fluids1'); if ~any(Iq); handles.fluids1=zeros(size(handles.Time)); 
handles.column{end+1}='fluids1'; end


I=[];
if ~any(strcmp(handles.column,'Mu1')); msgbox('Reduce the file first!'); return; end 

I(9)=find(strcmp(handles.column,'vel'),1,'last');
I(3)=find(strcmp(handles.column,'Mu1'),1,'last');
I(10)=find(strcmp(handles.column,'slip'),1,'last');
I(1)=find(strcmp(handles.column,'Time'),1,'last');
I(13)=find(strcmp(handles.column,'Rate'),1,'last');
I(4)=find(strcmp(handles.column,'shear1'),1,'last');
I(2)=find(strcmp(handles.column,'Normal'),1,'last');
I(5)=find(strcmp(handles.column,'SlipVel_Enc_1'),1,'last');
I(6)=find(strcmp(handles.column,'Slip_Enc_1'),1,'last');
I(7)=find(strcmp(handles.column,'SlipVel_Enc_2'),1,'last');
I(8)=find(strcmp(handles.column,'Slip_Enc_2'),1,'last');
I(11)=find(strcmp(handles.column,'LVDT_1'),1,'last');
I(12)=find(strcmp(handles.column,'LVDT_2'),1,'last');
I(14)=find(strcmp(handles.column,'fluids1'),1,'last');

% prova preallocazione
C=zeros(100000,1);
N=zeros(100000,1);
M=zeros(100000,1);
S=zeros(100000,1);

    for j=1:length(I); %1:length(handles.column)
        C(j,1)={'%10.6f '};
        if j==length(I)
            C(j,1)={'%10.6f\n'};
        end
    end
    C1=cell2mat(C');
    
    for j=1:length(I); %1:length(handles.column) 
        N(j,1)={['handles.' handles.column{I(j)} '(l,1),']};
        if j==length(I) 
             N(j,1)={['handles.' handles.column{I(j)} '(l,1)']};
        end
    end
    N1=cell2mat(N');
    
    
    for j=1:length(I); %1:length(handles.column) 
        M(j,1)={['''' handles.column{I(j)} '''' ',']};
        if j==length(I) 
             M(j,1)={['''' handles.column{I(j)} '''']};
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
if ~ strcmp(fieldnames(handles),'dt'); msgbox('ATTENTION: handles.dt=none'); end

%% OPEN (apertura file .txt con header)
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd('C:\LHVRFA\Data\Stefano')
cd txt
%ripulisci precedente

ah_=get(handles.axes1,'children'); set(ah_,'XData',[],'YData',[]);
ah_=get(handles.axes2,'children'); set(ah_,'XData',[],'YData',[]);
ah_=get(handles.axes3,'children'); set(ah_,'XData',[],'YData',[]);

I=strcmp(fieldnames(handles),'column');

if any(I) % any: Determine whether any array elements are nonzero
    for i=1:length(handles.column)
eval(['handles=rmfield(handles,''', handles.column{i}, ''');']) % rmfield: Remove fields from structure
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

[FileName,PathName] = uigetfile('*.txt*','C:\LHVRFA\Data\Stefano\txt');

file = [FileName,PathName];
cd (PathName)

%definisce i parametri da matrice
%handles.column=importdata(FileName,'\t',1);

fid=fopen(FileName,'r');
for i=1:3
file1=fgets(fid); % fgets: Read line from file, keeping newline characters
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
        if any(strcmp(fieldnames(handles),'column')) & ...
                strcmp(A,2); handles.column{i-1}={[char(handles.column(i-1)), '2']}; 
        else handles.column(i)={A};
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
            if ~ischar(tline), break, end
            [I]=find(tline==char(44));
            if ~isempty(I); tline(I)=char(46); end
            file1.data(i,:)=sscanf(tline,'%f');
    end 
    else
        file1=importdata(FileName,'\t',3);
    end
fclose(fid);

% %% prova da csv
% 
% f1 = fopen(file,'r');  % Open text file
% param1=textscan(f1, '%*s %f %f %f %f %f %f %f',1,'Delimiter',',','HeaderLines',11);
% fclose(f1);
% f2 = fopen(file,'r');  % Open text file
% param2=textscan(f2, '%*s %f %f %f %f %f %f %f',1,'Delimiter',',','HeaderLines',12);
% fclose(f2);
% f3 = fopen(file,'r');  % Open text file
% data=textscan(f3, '%f %f %f %f %f %f %f %f', 'Delimiter',',','HeaderLines',14); %read the file considering header the first 14 rows
% fclose(f3);
% f4 = fopen(file,'r');  % Open text file
% freq=textscan(f4, '%*s %f',1, 'Delimiter',',','HeaderLines',5);
% fclose(f4);
% 
% handles.Timestamp = repmat(1000/freq{1,1},size(data{1,1},1),1); %crea il vettore dalla freq campionamento (ms)
% handles.Axial = data{1,2}; %spinta ass
% handles.LVDT = data{1,3}; %lvdt
% handles.Torque = data{1,4}; %torque
% handles.Angle = data{1,5}; %angle
% handles.Encoder2 = data{1,6}; %tachometer
% handles.Encoder = data{1,7}; %puls counter
% handles.FG = data{1,8}; %function generator
% 
% %% fine prova

h_=findobj('Tag','dt_value');

[ndt,vdt]=grp2idx(file1.data(:,1));
if numel(vdt) > 1; handles.dt=str2double(vdt(2)); 
else
    handles.dt=str2double(vdt(1))
end
    

%set(h_,'String',handles.dt);

handles.filename=FileName;

handles.sm=0;
handles.triggered=0;
handles.cutted=[0 0];
handles.loadT=0;
handles.shearT=0;
ll=1;
nn=length(file1.data(:,1));

%primo step:togliere tutto quello che ha un campionamento diverso da dt
%handles.rate=0:handles.dt:(length(file1.data)-1)*handles.dt;
%memorizza anche gli originali
%eval(['handles.' handles.column{1} ' = cumsum(file1.data(ll:nn,1)); '])
%eval(['handles.v' num2str(1) ' = cumsum(file1.data(ll:nn,1)); '])
handles.column{1}='Time';
num=length(handles.column);

for n=2:num
test=double(handles.column{n});   
if any(test==32); I=find(test~=32); handles.column{n}=char(test(I)); end
eval(['handles.' handles.column{n} ' = file1.data(ll:nn,' num2str(n) ');'])
end


handles.column{num+1}='Stamp';
eval(['handles.' handles.column{num+1} '= file1.data(ll:nn,1); '])

num=length(handles.column);
handles.column{num+1}='Rate';
eval(['handles.' handles.column{num+1} '= [1:1:length(file1.data(ll:nn,1))]''; '])


num=length(handles.column);
handles.column{num+1}='RateZero';
eval(['handles.' handles.column{num+1} '= [1:1:length(file1.data(ll:nn,1))]''; '])


handles.TimeZero=cumsum(handles.Stamp);
handles.Time=zeros(size(handles.Stamp));
handles.Time(1)=handles.Rate(1)*handles.Stamp(1);
handles.Time(2:end)=handles.Rate(1)*handles.Stamp(1) +cumsum(handles.Stamp(2:end)); %plotto il numero di riga


guidata(hObject, handles);
handles.zoom=0;

plotta_ora(handles);

cd ../
%% --- Funzione plotta_ora (exec on calling)

function plotta_ora(handles)

h_=findobj('Tag','Rate');
k_=findobj('Tag','slipON');

stato= get(h_,'Value') ;
statok=get(k_,'Value') ;

if stato==1 && statok==0; set(h_,'String','time(s)'); handles.X=handles.Time/1000;
elseif stato==1 && statok==1; set(h_,'String','slip(m)'); handles.X=handles.slip;
else
set(h_,'String','rate'); handles.X= handles.Rate;  %se non voglio anche il rate triggerato
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


%% PRINT
function PrintMenuItem_Callback(~, ~, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

%% CLOSE
function CloseMenuItem_Callback(~, ~, handles)
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


%% --- Executes on button press in zoom.

function zoom_Callback(hObject, ~, handles)
% hObject    handle to zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pippo=[handles.axes1, handles.axes2,handles.axes3];
%zoom xon;
h_=findobj('Tag','Rate');
stato=get(h_,'Value');
if stato==1
t_cut=handles.Time/1000;
else
t_cut=handles.Rate;
end

linkaxes(pippo,'x');
[xi,yi]=ginput(2) ;

mat(1,:)=abs(t_cut-ones(size(t_cut))*xi(1));
mat(2,:)=abs(t_cut-ones(size(t_cut))*xi(2));
ll(:,1)=find(mat(1,:)==min(mat(1,:)));
ll(:,2)=find(mat(2,:)==min(mat(2,:)));

%for i=1:length(handles.column)
  
%eval(['handles.', handles.column{i}, '=handles.' ...
%    handles.column{i}, '(ll(1,1):ll(1,2),1);'])
%end

h_=findobj('Tag','Rate');
stato= get(h_,'Value') ;
if stato==1; set(h_,'String','time(s)'); handles.X=handles.Time/1000;
else
set(h_,'String','rate'); handles.X= handles.Rate;  %se non voglio anche il rate triggerato
end

set(handles.axes1,'Xlim',[handles.X(ll(1,1)) handles.X(ll(1,2))]);
set(handles.axes2,'Xlim',[handles.X(ll(1,1)) handles.X(ll(1,2))]);
set(handles.axes3,'Xlim',[handles.X(ll(1,1)) handles.X(ll(1,2))]);
handles.zoom=1;

guidata(hObject, handles);

plotta_ora(handles)
linkaxes(pippo,'off');


%% --- Executes on button press in zoom_back.
function zoom_back_Callback(hObject, ~, handles)
% hObject    handle to zoom_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h_=findobj('Tag','Rate');
stato= get(h_,'Value') ;
if stato==1
    set(h_,'String','time(s)')
handles.X=handles.Time/1000;
else
set(h_,'String','rate')
handles.X= handles.Rate;  %se non voglio anche il rate triggerato
end
set(handles.axes1,'XLim',[handles.X(1) handles.X(end)]);
set(handles.axes2,'XLim',[handles.X(1) handles.X(end)]);
set(handles.axes3,'XLim',[handles.X(1) handles.X(end)]);

handles.zoom=0;
plotta_ora(handles)
guidata(hObject, handles);


%% EDIT1 (menu di scelta variabile da plottare,1)
function edit1_Callback(hObject, ~, handles)
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

h_=findobj('Tag','Rate');
stato= get(h_,'Value') ;

if stato==1; set(h_,'String','timems)'); handles.X=handles.Time/1000;
else
set(h_,'String','rate'); handles.X= handles.Rate;  %se non voglio anche il rate triggerato
end



eval(['plot(handles.X,handles.' handles.column{(handles.g1)} ',''ob'',''parent'',handles.axes1);']); 
legend(handles.axes1,[handles.column{handles.g1}])

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(~, ~, ~)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.


%% EDIT2
function edit2_Callback(hObject, ~, handles)
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

h_=findobj('Tag','Rate');
stato= get(h_,'Value') ;

if stato==1; set(h_,'String','time(s)'); handles.X=handles.Time/1000;
else
set(h_,'String','rate'); handles.X= handles.Rate;  %se non voglio anche il rate triggerato
end


eval(['plot(handles.X,handles.' handles.column{(handles.g2)} ',''ob'',''parent'',handles.axes2);']); 
legend(handles.axes2,[handles.column{handles.g2}])
% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, ~, ~)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% EDIT3
function edit3_Callback(hObject, ~, handles)
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


h_=findobj('Tag','Rate');
stato= get(h_,'Value') ;

if stato==1; set(h_,'String','time(s)'); handles.X=handles.Time/1000;
else
set(h_,'String','rate'); handles.X= handles.Rate;  %se non voglio anche il rate triggerato
end

eval(['plot(handles.X,handles.' handles.column{(handles.g3)} ',''ob'',''parent'',handles.axes3);']); 
legend(handles.axes3,[handles.column{handles.g3}])

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, ~, ~)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




%% --- Executes on button press in Rate.
function Rate_Callback(hObject, ~, handles)
% hObject    handle to Rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% hObject    handle to Rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Rate
guidata(hObject, handles);

plotta_ora(handles)


%% --- Executes on button press in offset.
function offset_Callback(hObject, ~, handles)
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

h_=findobj('Tag','Rate');
stato=get(h_,'Value');
if stato==1
t_cut=handles.Time/1000;
else
t_cut=handles.Rate;
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




%% --- Executes on button press in trigger.
function trigger_Callback(hObject, ~, handles)
% hObject    handle to trigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h_=findobj('Tag','Rate');
stato=get(h_,'Value');
if stato==1
t_cut=handles.Time/1000;
else
t_cut=handles.Rate;
end


A=find(strcmp(fieldnames(handles),'shearT'));
if isempty(A) 
    [xi,yi]=ginput(1) ;
mat(1,:)=abs(t_cut-ones(size(t_cut))*xi(1));
ll(:,1)=find(mat(1,:)==min(mat(1,:)));
else
prev_trig=find(handles.RateZero==handles.shearT);
k=menu(['trigger is' num2str(handles.shearT) '. Is that ok?'],'si','no');

if (k==1)
    ll(:,1)=prev_trig;
elseif (k==2)
[xi,yi]=ginput(1) ;
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

h_=findobj('Tag','Rate');
set(h_,'Value',1)


ax_=get(handles.axes1,'Children'); dataY=get(ax_,'YData'); dataX=get(ax_,'XData');
set(ax_,'XData',handles.Time,'YData',dataY); set(handles.axes1,'XLim',[handles.Time(1) handles.Time(end)]); 
ax_=get(handles.axes2,'Children'); dataY=get(ax_,'YData'); dataX=get(ax_,'XData');
set(ax_,'XData',handles.Time,'YData',dataY); set(handles.axes1,'XLim',[handles.Time(1) handles.Time(end)]); 
ax_=get(handles.axes3,'Children'); dataY=get(ax_,'YData'); dataX=get(ax_,'XData');
set(ax_,'XData',handles.Time,'YData',dataY); set(handles.axes1,'XLim',[handles.Time(1) handles.Time(end)]); 

guidata(hObject, handles);
plotta_ora(handles)




%% --- Executes on button press in decimate.
function decimate_Callback(hObject, ~, handles)
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

%% --- Executes on button press in cut_dt

function cut_dt_Callback(hObject, eventdata, handles)
% hObject    handle to cut_dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cut_dt as text
%        str2double(get(hObject,'String')) returns contents of cut_dt as a double
set(hObject,'Enable','on')

h_=findobj('Tag','Rate');
stato=get(h_,'Value');
if stato==1
t_cut=handles.Time/1000;
else
t_cut=handles.Rate;
end


handles.dt=str2double(get(hObject,'String'));
%set(hObject,'String',handles.dt,'BackgroundColor',[0.75 0.75 0.75]) 
    ll=find(handles.Stamp(:,1)==handles.dt); %,1,'first');
    if length(ll) <= 100; h=msgbox('attention: number of residuals less than 100'); waitfor(h); return; end
    %nn=find(handles.Stamp(:,1)==handles.dt,1,'last');

    for n=1:length(handles.column)
    eval(['handles.' handles.column{n} ' = handles.' handles.column{n} '(ll,1);'])
    end


ax_=get(handles.axes1,'Children'); dataY=get(ax_,'YData'); dataX=get(ax_,'XData');
set(ax_,'XData',dataX(ll),'YData',dataY(ll)); set(handles.axes1,'XLim',[dataX(ll(1)) dataX(ll(end))]); 
ax_=get(handles.axes2,'Children'); dataY=get(ax_,'YData'); dataX=get(ax_,'XData');
set(ax_,'XData',t_cut(ll),'YData',dataY(ll)); set(handles.axes2,'XLim',[t_cut(ll(1)) t_cut(ll(end))]); 
ax_=get(handles.axes3,'Children'); dataY=get(ax_,'YData'); dataX=get(ax_,'XData');
set(ax_,'XData',t_cut(ll),'YData',dataY(ll)); set(handles.axes3,'XLim',[t_cut(ll(1)) t_cut(ll(end))]); 

   
set(hObject,'String','cut_dt')
guidata(hObject, handles);
plotta_ora(handles)

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



%% --- Executes on button press in cut.
function cut_Callback(hObject, eventdata, handles)
% hObject    handle to cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h_=findobj('Tag','Rate');
stato=get(h_,'Value');
if stato==1
t_cut=handles.Time/1000;
else
t_cut=handles.Rate;
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
eval(['handles.' handles.column{n} ' =handles.' handles.column{n} '(ll(:,1):ll(:,2),:);'])
end


guidata(hObject, handles);
plotta_ora(handles)



%% --- Executes on button press in fft.
function fft_Callback(hObject, eventdata, handles)
% hObject    handle to fft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h_=findobj('Tag','Rate');

stato=get(h_,'Value');
if stato==1
t_cut=handles.Time/1000;
else
t_cut=handles.Rate;
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
keyboard
set(hObject,'Value',0)

%eval(['plot(handles.Timestamp,handles.' handles.column{s} ...
%    ',''ob'',''parent'',handles.axes' num2str(finestra) ');']); 
%eval(['legend(handles.axes' num2str(finestra) ... 
%    ',[handles.column{' num2str(s) '}])'])



%% SAVE--------------------------------------------------------------------
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cd('C:\LHVRFA\Data\Stefano')

[nome,pat]=uiputfile( ...
      {'*.m;*.fig;*.mat;*.mdl', 'All MATLAB Files (*.m, *.fig, *.mat, *.mdl)'; ...
        '*.*',                   'All Files (*.*)'}, ...
         'Save as',['~/',handles.filename]);

cd(pat)

for j=1:length(handles.column)
        if j==length(handles.column) 
            M(j,1)={['''' handles.column{j} '''']};
        else
            M(j,1)={['''' handles.column{j} '''' ',']};
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
%eval(['save(nome3,''-struct'',''handles'',' O1 ');'])
%save('parametri','-struct','handles','loadT','shearT','triggered','cutted'
%,'dt','sm')


%% FILE --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




%% INTERACTIVE--------------------------------------------------------------------
function Interactive_Callback(hObject, eventdata, handles)
% hObject    handle to Interactive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
keyboard
return



%% --- Executes on button press in calibration.
function calibration_Callback(hObject, ~, handles)
% hObject    handle to calibration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%esclude gli smooth

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

%elimina vecchi new.qualcosa

if any(strcmp(fieldnames(handles),'new'));
    nomi=handles.new;
    for i=1:length(nomi)
        handles = rmfield(handles, nomi(i));
        K=find(~strcmp(handles.new(i),handles.column));
        handles.column=handles.column(K);
    end
    handles=rmfield(handles,'new');
end

%% cancella i vecchi vettori salvati (in handle.nome)
% for i=12:size(nomi,1)
%     handles=rmfield(handles, nomi{i,1});
% end
    
%input diametri

dint=findobj('Tag','Rint'); rint=str2double(get(dint,'String'))/2000;
dext=findobj('Tag','Rext'); rext=str2double(get(dext,'String'))/2000;

% dint=str2double(get(dint,'String'));
% dext=str2double(get(dext,'String'));

%% calibration coefficients and calculations

ang_=findobj('Tag','CheckAngle');
enc_=findobj('Tag','CheckEncoder');

stato_ang_=get(ang_,'Value') ;
stato_enc_=get(enc_,'Value') ;

%corrections old datasets
contents=get(hObject,'Value');

%cal.e2 is an average between the two corrections for internal and external
%velocity input

if contents==2; cal.e2=0.998; cal.e=1.49997; end
if contents==3; cal.e2=0.998; cal.e=0.99997; end
if contents==4; cal.e2=1; cal.e=1.030897; end
if contents==5; cal.e2=1; cal.e=1; end
% calculations

cal.torque(1:10)=3/2/pi/(rext^3-rint^3)*1E-6; % è ok

cal.ax(1:10)=(100*9.807)/(pi*(rext^2-rint^2))*1E-6; % è ok

cal.lv(1:10)=-0.001; % è ok

cal.enc(1:10)=cal.e*(4/3*pi*((rext^2+rint*rext+rint^2)/(rext+rint))); % WIP

cal.enc2(1:10)=cal.e2*(4/3*pi*((rext^2+rint*rext+rint^2)/(rext+rint))/60); % è ok

cal.ang(1:10)=(4/3*pi*((rext^2+rint*rext+rint^2)/(rext+rint)))*1/360; %deg to no of turns

%% ELABORAZIONE DATI

% Normal stress
b=(strfind(handles.column,'Axial')); j=0; n=[];
for i=1:length(b); if ~isempty(b{i}); j=j+1; n(j)=i; end; end
for j=1:length(n); 
eval(['new.Normal=handles.' handles.column{n(j)} '*cal.ax(j);'])
end

% Shear, Mu
b=(strfind(handles.column,'Torque')); j=0; n=[];
for i=1:length(b); if ~isempty(b{i}); j=j+1; n(j)=i; end; end
for j=1:length(n); 
eval(['new.Shear=handles.' handles.column{n(j)} '*cal.torque(j);'])
eval(['new.Mu=new.Shear./new.Normal;']);
end

% Vertical displacement-thickness (if zeroed)
b=(strfind(handles.column,'LVDT')); j=0; n=[];
for i=1:length(b); if ~isempty(b{i}); j=j+1; n(j)=i; end; end
for j=1:length(n); 
eval(['new.v_disp=handles.' handles.column{n(j)} '*cal.lv(j);'])
end

% Slip: Encoder OR Angle depends on selection

if stato_ang_==1 && stato_enc_==0;
b=(strfind(handles.column,'Angle')); j=0; n=[];
for i=1:length(b); 
    if ~isempty(b{i}); j=j+1; n(j)=i; 
    end; 
end
for j=1:length(n); 
eval(['new.Slip=handles.Angle*cal.ang(j);'])
end
elseif stato_ang_==0 && stato_enc_==1;
b=(strfind(handles.column,'Encoder')); j=0; n=[];
for i=1:length(b); 
    if ~isempty(b{i}); j=j+1; n(j)=i; 
    end; 
end
for j=1:length(n); 
eval(['new.Slip=handles.Encoder*cal.enc(j);'])    
end
end

if stato_enc_==1 && stato_ang_==0 %calculate true angle and veq_rpm only if there is encoder selected
% True angle
cal.angle(1:10)=0.0174532925; %deg to radians!
b=(strfind(handles.column,'Angle')); j=0; n=[];
for i=1:length(b); if ~isempty(b{i}); j=j+1; n(j)=i; end; end
for j=1:length(n); 
eval(['new.tr_angle=sin(handles.' handles.column{n(j)} '*cal.angle(j));'])
end

%Equivalent velocity from rpm
b=(strfind(handles.column,'Encoder2')); j=0; n=[];
for i=1:length(b); if ~isempty(b{i}); j=j+1; n(j)=i; end; end
for j=1:length(n); 
eval(['new.veq_rpm=handles.' handles.column{n(j)} '*cal.enc2(j);'])
end
end

% % % if stato_enc_==0 && stato_ang_==1
% % % 
% % % %Slip vel from Angle vs Time interpolation
% % % b=(strfind(handles.column,'Angle')); j=0; n=[];
% % % for i=1:length(b); if ~isempty(b{i}); j=j+1; n(j)=i; end; end
% % % for j=1:length(n); 
% % % clear d0 d1 d2 v vel
% % % eval(['d0=handles.' handles.column{n(j)} ';'])
% % % 
% % % 
% % % 
% % % %che i massimi siano fuori dal rumore
% % % J=local_max(d0); I=[];
% % % if length(J) > 10
% % % [c,d]=hist(diff(d0(J))); Jc=find(c==max(c));
% % % r=abs(diff(d0(J))-d(Jc(1))); In=find(r>mode(diff(d))); I=J(In);
% % % end
% % % 
% % % I=J;
% % % % figure(2); plot(d0); hold on; plot(I,d0(I),'*r')
% % % 
% % % 
% % % if ~isempty(I)
% % % for ii=1:length(I)
% % %     if ii==length(I) & I(ii) + 2 < length(d0)
% % % d0(I(ii)+1:end)=d0(I(ii)+1:end)-d0(I(ii)+2)+d0(I(ii));
% % %     elseif ii==length(I) & I(ii) +2 >= length(d0)
% % %         break
% % %     else
% % % d0(I(ii)+1:I(ii+1))=d0(I(ii)+1:I(ii+1))-d0(I(ii)+2)+d0(I(ii));
% % %     end
% % % end
% % % end
% % % 
% % % %I=find(d0<0);
% % % %d0(I)=0;
% % % 
% % % 
% % % h=handles.Stamp/1000;
% % % d1=d0*cal.ang(1);
% % % 
% % % 
% % % %primo tipo di operazione
% % % % tipo=input('flt type: v0 v1 v2 v3 v4 ','s')
% % % % tipo='v2';
% % % % [v,v4,d00,d4]=smuut_o(d1,h,tipo);
% % % % y=d1;
% % % 
% % % 
% % % %altra opzione, smoothing dello slip con minimizzazione defasamento
% % % % y=interp_Enc(handles.Time,d1,j);
% % % % bb=zeros(size(y));
% % % % bb(2:end)=diff(y); bb(1)=bb(2);
% % % % v=bb./h;
% % % 
% % % %terzo tipo, no filtro
% % % test=1;
% % % if test==1
% % % L=2*floor(length(d1)/100)+1; 
% % % sigma=3;
% % % if j==2; L=1500; sigma=50; end
% % % sm = fspecial('gaussian',[L 1],sigma); % gaussian kernel where s= size of contour
% % % y = conv(d1, sm); d_sm=y(floor(L/2):length(d1)+floor(L/2)-1);
% % % y=[]; y=d_sm;
% % % else
% % % d_sm=d1;
% % % y=d_sm;
% % % end
% % % bb=zeros(size(d_sm));
% % % bb(2:end)=diff(d_sm); bb(1)=bb(2);
% % % v=bb./h;
% % % %if j==2; v=smooth(v,300); end
% % % 
% % % 
% % % 
% % % if j==1; del=zeros(size(y)); end
% % % if j==2; del0=del; end
% % % del=(max(del,y)); 
% % % if j==2; adel=abs(y-del)./del.*100; Iadel=find(adel < 1); 
% % % del(Iadel)=del0(Iadel);
% % % end
% % % 
% % % if j==1; eval(['new.SlipVel_Enc_' num2str(j) '=v']);
% % % eval(['new.Slip_Enc_' num2str(j) '=(d1);']);
% % % end%'=(dpippo./dtime); ']')
% % % 
% % % 
% % % if j==2; eval(['new.SlipVel_Enc_' num2str(j) '=v']); 
% % % eval(['new.Slip_Enc_' num2str(j) '=(d_sm);']);
% % % end%'=(dpippo./dtime); ']')
% % % 
% % % 
% % % bb=zeros(size(del));
% % % bb(2:end)=diff(del); bb(1)=bb(2);
% % % 
% % % %if j==1; Ivmax=find(v<=0.05); new.vel=zeros(size(v)); new.vel(Ivmax)=v(Ivmax); end
% % % %if j==2; Ivmax=find(v>0.05);  new.vel(Ivmax)=new.SlipVel_Enc_2(Ivmax); end
% % % new.slip=del;
% % % new.vel=bb./h;
% % % 
% % % end 
% % % end
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


%% diametro interno 
function Rint_Callback(~, ~, ~)
% hObject    handle to Rint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Rint as text
%        str2double(get(hObject,'String')) %returns contents of Rint as a double

% --- Executes during object creation, after setting all properties.
function Rint_CreateFcn(hObject, ~, ~)
% hObject    handle to Rint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% diametro esterno
function Rext_Callback(~, ~, ~)
% hObject    handle to Rext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Rext as text
%        str2double(get(hObject,'String')) %returns contents of Rext as a double


% --- Executes during object creation, after setting all properties.
function Rext_CreateFcn(hObject, ~, ~)
% hObject    handle to Rext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




%% --- Executes on button press in Gefran.
% function Gefran_Callback(hObject, eventdata, handles)
% % hObject    handle to Gefran (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hint: get(hObject,'Value') returns toggle state of Gefran
% stato=get(hObject,'Value');
% if stato==1  
%     ButtonName = questdlg('GEFRAN preliminary operations: (1) identify the main sampling rate with cut_dt (2) decimate the data to 1s', ...
%                          'Do you want to proceed now?', 'No', 'Yes','No');
%                     
%    switch ButtonName,
%      case 'No',
%      disp('no');
%      case 'Yes',
%    
% I=find(~strcmp(handles.column,'SpeedGEF')); handles.column=handles.column(I);
% I=find(~strcmp(handles.column,'timeGEF')); handles.column=handles.column(I);
% I=find(~strcmp(handles.column,'TorqueGEF')); handles.column=handles.column(I);
% 
% Path2='/media/disk1/shivadir/Shiva Experiments';
% name=handles.filename(1:4);
% 
% list=dir([Path2 '/' name '*']);
%     if ~isempty(list); Path2=[Path2 '/' list.name '/']; 
%     name2=dir([Path2 '*.txt']);
%     FileGEF=name2.name;
%     else
%     [FileGEF,Path2] = uigetfile( '/media/disk1/shivadir/*.txt', ...
%     'Multiple File Detected: Select the file GEFRAN to load');
%     name2='';
%     end
% 
%     
% gefran1=1; k=0 ;   
% while ~isstruct(gefran1)
%     k=k+1;
%     gefran1=importdata([Path2 '/' FileGEF],'\t',k);
% end
% 
% I=find(strncmp(gefran1.textdata,'Time	Speed',6)); if ~isempty(I); new.timeGEF=gefran1.data(:,1); new.VGEF(:,1)=gefran1.data(:,2); end
% I=find(strncmp(gefran1.textdata,'Act Torque',6)); if ~isempty(I); new.TqGEF(:,1)=gefran1.data(:,I(1)-1); end
% I=find(strncmp(gefran1.textdata,'Speed',5)); if ~isempty(I); new.VGEF(:,1)=gefran1.data(:,I(1)-1); end
% I=find(strncmp(gefran1.textdata,'time',4)); if ~isempty(I); new.timeGEF(:,1)=gefran1.data(:,I(1)-1); end
% 
% dt=diff(new.timeGEF); dt(end+1)=dt(1);
% 
% %for j=1:length(handles.column);
% %    eval(['new.' handles.column{j} '= downsample(handles.' handles.column{j} ',' num2str(25) ');'])
% %end
% %dati di calibrazione
% 
% nomi=[];
% [a,b]=size(handles.column); [aa,bb]=size(fieldnames(new));
% inp1=handles.column;
% if aa==1 & aa==b | bb==1 & a==1; inp1=handles.column'; end
% nomi=[inp1 ; fieldnames(new)];
% 
% handles.column=[];
% handles.column=nomi';
% nomi2=fieldnames(new);
% 
% for k=1:length(nomi2)
%  eval(['handles.' char(nomi2(k)) '=new.' char(nomi2(k)) ';'])
% end
% 
% guidata(hObject, handles);
% 
%    end
% end % switch
% 
% %%%%%%%%%%% fine GEF


%% --- Executes on button press in running.
function running_Callback(hObject, ~, handles)
% hObject    handle to running (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h_=findobj('Tag','Rate');
stato=get(h_,'Value');
if stato==1
t_cut=handles.Time/1000;
else
t_cut=handles.Rate;
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
eval(['plot(handles.Rate,handles.' handles.column{s} ',''ob'',''parent'',handles.axes' num2str(finestra) ''');']); 
set(ax_,'Xlim',posx)

guidata(hObject, handles);
plotta_ora(handles)



%% --- Executes on button press in fluid.
function fluid_Callback(~, ~, ~)
% hObject    handle to fluid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fluid



%% PRINT --------------------------------------------------------------------
function print_Callback(~, ~, handles)
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


%ax_=get(handles.axes1,'Parent');
%ah_=get(ax_,'Children')
%hf=figure;
%for n=1:length(ah_)
%eval(['h_new = copyobj(ah_(' num2str(n) '),hf);']);
%end




%% --- Executes on button press in outlayers.
function outlayers_Callback(hObject, ~, handles)
% hObject    handle to outlayers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%zoom xon;
h_=findobj('Tag','Rate');
stato=get(h_,'Value');
if stato==1
t_cut=handles.Time/1000;
else
t_cut=handles.Rate;
end

button=1; i=0;
while button==1
i=i+1;
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




%% --- Executes on selection change in edit1LB.
function edit1LB_Callback(hObject, ~, handles)
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



% --- Executes during object creation, after setting all properties.
function edit1LB_CreateFcn(hObject, ~, ~)
% hObject    handle to edit1LB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




%% --- Executes on selection change in edit2LB.
function edit2LB_Callback(hObject, ~, handles)
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


% --- Executes during object creation, after setting all properties.
function edit2LB_CreateFcn(hObject, ~, ~)
% hObject    handle to edit2LB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Executes on selection change in edit3LB.
function edit3LB_Callback(hObject, ~, handles)
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



% --- Executes during object creation, after setting all properties.
function edit3LB_CreateFcn(hObject, ~, ~)
% hObject    handle to edit3LB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






%% FIGURE --------------------------------------------------------------------
function Figure_Callback(~, ~, ~)
% hObject    handle to Figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% LOAD --------------------------------------------------------------------
function Load_Callback(hObject, ~, handles)
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




[FileName,PathName] = uigetfile('C:\Users\Aretusini\Dropbox\Dottorato\ROSA');

cd (PathName)

data=load(FileName);


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

handles.column=fieldnames(data);

for i=1:length(handles.column)
    eval(['handles.' handles.column{i} '=data.' handles.column{i} ';'])
end


h_=findobj('Tag','edit1LB'); set(h_,'String',handles.column);
h_=findobj('Tag','edit2LB'); set(h_,'String',handles.column);
h_=findobj('Tag','edit3LB'); set(h_,'String',handles.column);


guidata(hObject, handles);
handles.zoom=0;

plotta_ora(handles);






%% BINARY
function binary_Callback(~, ~, handles)
% hObject    handle to binary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[nome,pat]=uiputfile( ...
      {'*.txt', 'All MATLAB Files (*txt)'; ...
        '*.*',                   'All Files (*.*)'}, ...
         'Write as',['~/',handles.filename]);
cd (pat)
C=zeros(10000,1);
N=zeros(10000,1);
M=zeros(10000,1);
S=zeros(10000,1);

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



% --- Executes on button press in slipON.
function slipON_Callback(hObject, eventdata, handles)
% hObject    handle to slipON (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of slipON




% --- Executes on button press in CheckAngle.
function CheckAngle_Callback(hObject, eventdata, handles)
% hObject    handle to CheckAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckAngle


% --- Executes on button press in CheckEncoder.
function CheckEncoder_Callback(hObject, eventdata, handles)
% hObject    handle to CheckEncoder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckEncoder


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
