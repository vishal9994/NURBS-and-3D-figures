function NURBS(varargin)
if nargin == 0
  clear all
end

close all
clc

fprintf(['--------------------\n', ...
  'CAD Project under guidance of Dr. P.K. Jain Sir \n', ...
  '--------------------\n\n'])

fig = figure('Position', [120 100 1100 700],'NumberTitle','off','Name','CAGD : NURBS', 'color',[1 1 1]);
set(fig, 'MenuBar','none');
set(fig,'ToolBar','none');
hold on

set(fig,'WindowButtonMotionFcn',@move_mouse);
set(fig,'WindowButtonDownFcn',@click_mouse);
set(fig,'WindowScrollWheelFcn',@scroll_mouse);
set(fig,'KeyPressFcn',@keyboard )

global Deg KVect CPoint Weights s

% Default values
Deg = 3;
KVect = [0,1,1,1,2,3,4,5,6,6,6,7];%[0, 0, 0, 1/4, 1/4, 1/2, 1/2, 3/4, 3/4, 1, 1, 1]
s = 0;
CPoint = [];
Weights = [];
x_min = -10;
x_max = 10;
min_y = 0;
y_max = 10;
axis( [ x_min-1, x_max+1, min_y-1, y_max+1 ] )
%axis square
axis off

global Select point_move Labels Joint Index

Select = 1;
point_move = 0;
Labels = 1;
Joint = 0;
Index = [];


init_curve()

% -----------------------------------------------------------------------

function init_curve()
global Basis Deg KVect CPoint Weights Index

[~, Basis] = DeBoor(Deg,KVect,Deg+1,length(KVect)-Deg);
numCP = length(KVect) - Deg - 1;

if size(CPoint,1) < numCP
  AddCP = numCP - size(CPoint,1);
  disp(['Add ', num2str(AddCP), ' control points.'])
  set(gcf,'Pointer','CrossHair')
elseif size(CPoint,1) > numCP
  DelCP = size(CPoint,1) - numCP;
  CPoint = CPoint(1:end-DelCP,:);
  Index = [];
  Weights = Weights(1:end-DelCP);
  disp([':: Removed', num2str(DelCP), ' control points.'])
  update()
else
  update()
end

% -----------------------------------------------------------------------

function update()
global Basis Deg KVect CPoint Weights Index Labels Joint Select s
cla
% Control polygon
plot( CPoint(:,1), CPoint(:,2), '-.', 'Color', .4*[0 1 1], 'LineWidth', 1 );
% Control points
plot( CPoint(:,1), CPoint(:,2), 'bo', 'MarkerSize', 12, 'MarkerFaceColor', [0.000, 0.122, 0.247] );

switch Labels
  case 1
    % Control point labels
    for pp = 1:length(CPoint)
      text( CPoint(pp,1), CPoint(pp,2), ['    P_{', num2str(pp), '}'], 'FontSize', 12 )
    end
  case 2
    % Weights
    for pp = 1:length(CPoint)
      text( CPoint(pp,1), CPoint(pp,2), ['    W (', num2str(Weights(pp)), ')'], 'FontSize', 12 )
    end
end

denom = Weights* Basis;
x = (Weights.*CPoint(:,1)' * Basis) ./ denom;
y = (Weights.*CPoint(:,2)' * Basis) ./ denom;

if s==0
    for i = 1:length(x)
    plot(x(1:i),y(1:i),'Color', [44 133 72]/255, 'LineWidth', 3)
    drawnow
    end
    s = s+1;
else
plot( x, y, 'Color', [44 133 72]/255, 'LineWidth', 3 );
end
if Joint    
  Knots = unique(KVect((Deg+1):(end-Deg)));
  for xx = Knots
    Joint = For_t(xx);
    plot( Joint(1), Joint(2), 's', 'Color', [44 133 72]/255, 'MarkerSize', 10, 'MarkerFaceColor', [44 133 72]/255 )
    text( Joint(1), Joint(2)-.3, ['t = ', num2str(xx)], 'FontSize', 12 )
  end
end

if Select == 1
  plot( CPoint(Index,1), CPoint(Index,2), 'o', 'MarkerSize', 16, 'Color', .7*[1 1 1], 'LineWidth', 3 );
else
  plot( [CPoint(Index,1) CPoint(Index+1,1)], [CPoint(Index,2) CPoint(Index+1,2)], 'Color', .7*[1 1 1], 'LineWidth', 3 );
end

% -----------------------------------------------------------------------

function move_mouse(~,~)
global point_move move_along CPoint Index

if point_move
  CurrentPoint = get(gca,'CurrentPoint');
  point1 = CurrentPoint(1,1:2);  
 switch move_along
 case 0
 CPoint(Index,:) = point1;
 case 1
 CPoint(Index,1) = point1(1);
 case 2
 CPoint(Index,2) = point1(2);
  end
  
update()
end

% -----------------------------------------------------------------------

function click_mouse(~,~)
global Select point_move  CPoint Weights KVect Deg Index

% Left click
if strcmpi(get(gcf,'SelectionType'), 'Normal')
  
 if point_move
    point_move = 0;
    Index = [];
    disp(':: Releasing control point')
    update();
  else
    if strcmpi( get(gcf,'Pointer'),'CrossHair')
      
      CurrentPoint = get(gca,'CurrentPoint');
      Pt = CurrentPoint(1,1:2);      
      CPoint(end+1,:) = Pt;
      Weights(end+1) = 1;
      plot( CPoint(end,1), CPoint(end,2), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', [0.000, 0.122, 0.247] );
      disp(':: CP Added')
      
 if size(CPoint,1) == (length(KVect)-Deg-1)
    set(gcf,'Pointer','Cross')
    update()
 end
      
 end
end
  
  % Right click
elseif strcmpi(get(gcf,'SelectionType'), 'Alt')
  
if Select == 1
 disp('::CP Selected')
  Index = closest();
  update()
    
 elseif Select == 2
 disp(':: Edge Selected')
  Index = closest();
 update()
  end
  
end

% -----------------------------------------------------------------------

function scroll_mouse(~,press)
global Select Weights Index Deg KVect s

if Select == 1
  dd = 10;
  value = press.VerticalScrollCount / dd;
  Weights(Index) = Weights(Index) + value;
  if Weights(Index) < 0
      
    warning('Negative weight value!')
  end
update()
elseif Select == 2
  if (mod(Deg,2) == 0 && Select == 1 && length(Index)) || (mod(Deg,2) == 1 && Select == 2 && length(Index))
    dd = 10;
    value = press.VerticalScrollCount / dd;
    KVect = [KVect(1:round(Deg/2)+Index) KVect(round(Deg/2)+1+Index:end) + value*0.5];
    init_curve()
  end 
end

% -----------------------------------------------------------------------

function keyboard(~,press)
global KVect Deg Labels Joint CPoint Weights Index Select point_move move_along

switch press.Character
      
    case {'b', 'B'}
    % Plot the basis
    plot_basis();
    case {'d', 'D'}
    % Change the degree
    dd = input('Deg (e.g. quadratic = 2): ');
    if length(KVect) < (2 + 2*Deg)
      warning('Knot vector does not satisfy minimal length for this degree!')
    else
      Deg = dd;
      init_curve()
    end   
    case {'e', 'E'}
    % select mode to edges
    Select = 2;
    Index = [];
    disp(':: Select Mode [Edge]')
    update();    
    case {'m', 'M'}
    % move control point
    if (Select == 1 && length(Index))
      disp(':: Moving control point')
      point_move = 1;
      move_along = 0;
    end   
    case {'i', 'I'}
    % Display information about the curve
    kk = strtrim(sprintf('%.4g  ', KVect));
    ww = strtrim(sprintf('%.4g  ', Weights));
    fprintf(['\nKnotvector: [', kk, ']\n', ...
      'Deg: ', num2str(Deg), '\n', ...
      'Number of control points: ', num2str(size(CPoint,1)), '\n', ...
      'Weights: ', ww, '\n\n'])  
    case {'j', 'J'}
    % visibility of the knots
    Joint = ~Joint;
    update()
    case {'k', 'K'}
    % Change the knot vector
    KVect = input('knot (between [ and ]): ');
    init_curve()
    case {'l', 'L'}
    %0 no label,1 cp label,2 weight label)
    Labels = mod(Labels+1,3);
    switch Labels
      case 1
        disp('::control point labels')
      case 2
        disp('::Weights')
    end
    update()    
    case {'p', 'P'}
    % position of selected CP
    if Select == 1
      CPoint(Index,:) = input('Position (between [ and ]): ');
      update()
    end  
    case {'t', 'T'}
    % curve at given 't'
    Val = input('Evaluate curve at t: ');
    if Val < KVect(Deg+1) || Val >= KVect(end-Deg)
      warning('t value out of range!')
    end
    pp = For_t(Val);
    plot(pp(1), pp(2), 'rx', 'MarkerSize', 12, 'LineWidth', 3)
    text(pp(1), pp(2)-.3, ['t = ', num2str(Val)], 'Color', [1 0 0], 'FontSize', 12 )  
    case {'v', 'V'}
    % Select to Vertices
    Select = 1;
    Index = [];
    set(gcf,'Name','NURBS [Vertex]');
    disp(':: Select Mode [Vertex]')
    update(); 
    case {'w', 'W'}
    %enter weight of selected CP
    if Select == 1
      Weights(Index) = input('Weight: ');
      update();
    end
    case {'q', 'Q'}
    % Quit
    close all
    clear all
    clc
    case {'s', 'S'}
    % Save the curve in SVG format
    % http://www.mathworks.nl/matlabcentral/fileexchange/7401-scalable-vector-graphics-svg-export-of-figures
    plot2svg
    case {'h', 'H'}
    % Plot first derivative of the curve
    fderiv()
    case {'z', 'Z'}
    % Extrude
    %vector along which it has to be extruded
    vect=[];
    vect = input('vector[enter the vector]: ');
    extrude(vect)
    case {'a', 'A'}
    % Toggle axis visibility
    if strcmpi(get(gca, 'Visible'), 'On')
      axis off;
    else
      axis on;
    end
    case {'x', 'X'}
    %vector along which it has to be extruded
    info()
    case {'c', 'C'}
    %vector along which it has to be extruded
    surface1()
end

% -----------------------------------------------------------------------

function Pt = For_t(tvalue)
global KVect Deg CPoint Weights  
weight=Weights;
W = unique(Weights);
if(length(W))>=1
    start = Deg+1;
    endv = length(KVect)-Deg;
    p=weight'.*CPoint;
    p(:,3)=weight';
    Basis1=DeBoor1(Deg,KVect,start,endv,tvalue);
    Basis1=Basis1';
    cc=Basis1*p;
    cc= cc/cc(1,3);
    cc(:,3) = [];
    %CurrentPoints=S/S(1,3)
    cpp = cc;


  
end

Pt = cpp;

% -----------------------------------------------------------------------

function Index = closest()
global CPoint Select;

cpp = get(gca,'CurrentPoint');
Pt = cpp(1,1:2);

if Select == 1
  
  dist = zeros(1,length(CPoint));
  
  for P = 1:length(dist)
    dist(P) = norm( CPoint(P,:)-Pt, 2);
  end
  
else
  
  dist = zeros(1,length(CPoint)-1);
  ppx = Pt(1);
  ppy = Pt(2);
  
  for jj = 1:length(dist)
    aax = CPoint(jj,1);
    aay = CPoint(jj,2);
    bbx = CPoint(jj+1,1);
    bby = CPoint(jj+1,2);
    % Difference
    diffx = bbx - aax;
    diffy = bby - aay;
    Val = ((ppx-aax)*diffx + (ppy-aay)*diffy)/(diffx^2 + diffy^2);
    
    if Val < 0
      Val = 0;
    elseif Val > 1
      Val = 1;
    end
    
    dist(jj) = norm( (1-Val)*CPoint(jj,:)+Val*CPoint(jj+1,:)-Pt, 2);
  end
  
end

[~,Index] = min( dist );

% -----------------------------------------------------------------------

function fderiv()
global CPoint KVect Deg

npts = zeros(length(CPoint)-1,2);

for Q = 1:length(npts)
  npts(Q,:) = (Deg/(KVect(Q+1+Deg)-KVect(Q+1)))*(CPoint(Q+1,:)-CPoint(Q,:));
end

H = figure('NumberTitle','off','Name','First Derivative','color',[1 1 1]);
set(H, 'MenuBar','none');
set(H,'ToolBar','figure');
hold on
% Control polygon
plot( npts(:,1), npts(:,2), '-.', 'Color', .4*[0 1 1], 'LineWidth', 1 );
% Control points
plot( npts(:,1), npts(:,2), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', [0.000, 0.122, 0.247] );

[~, nbasis] = DeBoor(Deg-1,KVect(2:end-1),Deg,length(KVect)-Deg-1);

X = npts(:,1)' * nbasis;
Y = npts(:,2)' * nbasis;
X(:,end)=[];
Y(:,end)=[];
for i = 1:length(X)
    plot(X(1:i),Y(1:i),'Color', [230, 133, 14]/255, 'LineWidth', 3)
    drawnow
end

for Q = 1:length(npts)
  text( npts(Q,1), npts(Q,2), ['   Q_{', num2str(Q), '}'], 'FontSize', 12 )
end

axis off

% -----------------------------------------------------------------------

function plot_basis()
global KVect Deg CPoint Weights
if strcmpi(get(gcf,'CurrentModifier'), 'Shift')
  From = 1;
  To = length(KVect);
else
  From = Deg+1;
  To = length(KVect)-Deg;
end

[mat_t, basis] = DeBoor(Deg,KVect,From,To);
colors = ['8b327a'; '633131'; 'c70202'; 'ecbd00'; '3f9a02'; '1e927b'];

for ii = 1:length(CPoint)
  k = mod(ii,length(colors))+1;
  colour = hex2dec([colors(k,1:2); colors(k,3:4); colors(k,5:6)])/255;
  plot(CPoint(ii,1), CPoint(ii,2), 'o', 'MarkerSize', 16, 'Color', colour, 'LineWidth', 3 );
end
%Basis dimensions
breadth = 1100;
height = breadth / (length(CPoint) - Deg + 1) + 100;

fig = figure('NumberTitle','off','Name','Basis','color',[1 1 1]);
set(fig, 'MenuBar','none');
set(fig,'ToolBar','none');
set(fig,'Position',[125 100 breadth height]);
hold on
denom = Weights * basis;

for ii = 1:length(CPoint)
  k = mod(ii,length(colors))+1;
  colour = hex2dec([colors(k,1:2); colors(k,3:4); colors(k,5:6)])'/255;
  plot(mat_t, Weights(ii)*basis(ii,:) ./ denom, 'Color', colour, 'LineWidth', 2)
end
axis tight

% -----------------------------------------------------------------------

function extrude(vect)
global CPoint Weights Deg KVect
  From = Deg+1;
  To = length(KVect)-Deg;
[~, basis] = DeBoor(Deg,KVect,From,To);
x_value = vect(1,1);
y_value = vect(1,2);
z_value = vect(1,3);

denom = Weights* basis;
x = (Weights.*CPoint(:,1)' * basis) ./ denom;
y = (Weights.*CPoint(:,2)' * basis) ./ denom;
pp =(length(x));
cp =(length(CPoint));
%disp(x)
x_valuep=linspace(0,x_value,pp);
y_valuep=linspace(0,y_value,pp);
z_valuep=linspace(0,z_value,pp);

mat = ones(pp,4);
mat(:,1) = x;
mat(:,2) = y;
mat(:,3) = 0;

cmat = ones(cp,4);
cmat(:,1) = CPoint(:,1);
cmat(:,2) = CPoint(:,2);
cmat(:,3) = 0;

ff = linspace(0,z_value,pp);
fig = figure('NumberTitle','off','Name','Extrusion','color',[1 1 1]);
set(fig, 'MenuBar','figure');
set(fig,'ToolBar','figure');
set(fig,'Position',[120 100 700 500]);
hold on

% Control polygon
    plot3( CPoint(:,1), CPoint(:,2),cmat(:,3) ,'-.', 'Color', [1 0.2 0.4], 'LineWidth', 1 );
% Control points
    plot3( CPoint(:,1), CPoint(:,2),cmat(:,3), 'bo', 'MarkerSize', 3, 'MarkerFaceColor', [0.000, 0.122, 0.247] );
    
    tc_mat=[1,0,0,0;0,1,0,0;0,0,1,0;x_value,y_value,z_value,1];
    nc_mat=cmat*tc_mat;
    
 % Control points
    plot3( nc_mat(:,1),nc_mat(:,2) ,nc_mat(:,3), 'bo', 'MarkerSize', 3, 'MarkerFaceColor', [0.000, 0.122, 0.247] );
 % Control polygon
    plot3( nc_mat(:,1),nc_mat(:,2) ,nc_mat(:,3),'-.', 'Color', [1 0.2 0.4], 'LineWidth', 1 );
 view(-160,-40)
 
for E=1:length(ff)
    
    t_mat=[1,0,0,0;0,1,0,0;0,0,1,0;x_valuep(1,E),y_valuep(1,E),z_valuep(1,E),1];
    nm_mat=mat*t_mat;
    x = nm_mat(:,1);
    y = nm_mat(:,2);
    z = nm_mat(:,3);
    plot3( x,y,z,'Color', .4*[0 1 1], 'LineWidth',1 );
    drawnow
    view(-160,-40)
end

axis on
axis tight

% -----------------------------------------------------------------------

function surface1()
global CPoint Weights Deg KVect
  From = Deg+1;
  To = length(KVect)-Deg;
[~, basis] = DeBoor(Deg,KVect,From,To);
denom = Weights* basis;
x = (Weights.*CPoint(:,1)' * basis) ./ denom;
y = (Weights.*CPoint(:,2)' * basis) ./ denom;
%x(:,1) = [];
%x(:,end) = [];
%y(:,1) = [];
%y(:,end) = [];
pp =(length(x));
cp =(length(CPoint));
fig = figure('NumberTitle','off','Name','Revolve','color',[1 1 1]);
set(fig, 'MenuBar','figure');
set(fig,'ToolBar','figure');
set(fig,'Position',[120 100 700 500]);
hold on
mat = ones(pp,4);
mat(:,1) = x;
mat(:,2) = y;
mat(:,3) = 0;

cmat = ones(cp,4);
cmat(:,1) = CPoint(:,1);
cmat(:,2) = CPoint(:,2);
cmat(:,3) = 0;

jj=linspace(0,2*pi,pp-20);
mx=[];
my=[];
mz=[];
for p=1:length(jj)
    X = mat(:,1);
    Y = mat(:,2)*cos(jj)-mat(:,3)*sin(jj);
    Z = mat(:,2)*sin(jj)+mat(:,3)*cos(jj);
    
    %Cx = cmat(:,1);
    %Cy = cmat(:,2)*cos(jj)-cmat(:,3)*sin(jj);
    %Cz = cmat(:,2)*sin(jj)+cmat(:,3)*cos(jj);   
 
    %mx = cat(1,mx,X);
    %my = cat(1,my,Y);
    %mz = cat(1,mz,Z);
    plot3( X, Y, Z,'Color',[153 0 0]/255, 'LineWidth',1);
    %plot3( X, Y, Z, 'o', 'MarkerSize', 2, 'MarkerFaceColor', [0.000, 0.122, 0.247] );
    % Control points
    %plot3(Cx,Cy ,Cz, 'bo', 'MarkerSize', 3, 'MarkerFaceColor', [0.000, 0.122, 0.247] );
    % Control polygon
    %plot3( Cx,Cy ,Cz,'-.', 'Color', [1 0.2 0.4], 'LineWidth', 1 );
    
    
    view(-25,-28)
end
%mx(:,2)=my;
%mx(:,3)=mz;
%size(my)%mesh(mx);
kk=x;
radius = abs(y);
for q=1:length(kk)
    radii=radius(q);
    center = [kk(q),0,0];
    normal=[1,0,0];
    theta=0:0.01:2*pi;
    v=null(normal);
    points=repmat(center',1,size(theta,2))+radii*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
    plot3(points(1,:),points(2,:),points(3,:),'Color',[204 102 2]/255,'LineWidth',1);
    %plot3(points(1,:),points(2,:),points(3,:), 'bo', 'MarkerSize', 4, 'MarkerFaceColor', [0.000, 0.122, 0.247] );
    view(-25,-28)
    if cp<=5
    drawnow
    end
end
view(-25,-28)
axis on
axis tight

% -----------------------------------------------------------------------

function info()
Keyboard = {'A';'B';'C';'D';'E';'H';'I';'J';'K';'L';'M';'P';'Q';'S';'T';'V';'W';'Z'};
Application = {'Show Axis';'Plot Basis';'Revolve Curve';'Change Degree';'Select edges';
    'Plot First Derivative';'Curve Information';'Show Joints';'Change Knots Vector';'Show labels';'Move point';
    'Assign position to CP';'Quit';'Save plot as svg';'Evaluate curve at "t"';'Select vertices';
    'Assign weigths';'Extrude'};
T = table(Application,'RowNames',Keyboard);   
TString = evalc('disp(T)');
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');

FixedWidth = get(0,'FixedWidthFontName');

fig = figure('NumberTitle','off','Name','Instructions','color',[1 1 1]);
set(fig, 'MenuBar','none');
set(fig,'ToolBar','none');
set(fig,'Position',[1220 450 300 360]);
axis off
hold on
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
    'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);

% -----------------------------------------------------------------------
