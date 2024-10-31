%% level set methods
function LevelCal(method)
%-----------------------------------------------------------------------
nodesNum = 100;       % number of space gridpoints
dt = 8e-3;    % time step
tf = 0.6;    % final time
reinitNum = 100;       % number of reinitialization steps
nsteps = 10;  % number of steps with graphic output
%-----------------------------------------------------------------------
if nargin<1  % if the user didnt select the velocity type
    method = 2; 
end
nTimeStep = ceil(tf/dt); dt = tf/nTimeStep;
x = linspace(-1,1,nodesNum); y = linspace(-1,1,nodesNum); h = x(2)-x(1);
[X,Y] = meshgrid(x); ax = [min(x) max(x) min(y) max(y)];
%-----------------------------------------------------------------------
% initial phi 
P = zeros(nodesNum,nodesNum);
for i = 1:nodesNum
    for j = 1:nodesNum
        x(i) = -1 + (i-1) * h; 
        y(j) = -1 + (j-1) * h; 
        if isInside(x(i),y(j),0,0,1)
            P(i,j) = -1; 
        else 
            P(i,j) = 1; 
        end 
    end
end
%surf(x,y,P); hold on; axis ([ax]); 
%-----------------------------------------------------------------------
for it = 1:nTimeStep
    switch method
        case 1, F = 0.05;              % movement with constant velocity
        case 2, F = -curvature(P,h)*5e-3;    % movement under curvature
    end
        P = P-dt*FabsgradP(P,h,F);                        % level set update
    for ir = 1:reinitNum                               % reinitialization steps
        P = P-dt*FabsgradP(P,h,P./sqrt(P.^2+(2*h)^2),1);
    end

    if it==1 | floor(nsteps*it/nTimeStep)>floor(nsteps*(it-1)/nTimeStep) % visualization
        clf
        subplot(1,2,1), contourf(x,x,-P,[0 0],'k-') % this makes the graphs look better
        axis equal, axis(ax)
        title(sprintf('geometry at t=%0.2f',it*dt))
        subplot(1,2,2), surf(x,x,-P,'EdgeAlpha',.2), hold on
        patch([-1 1 1 -1],[-1 -1 1 1],[0 0 0 0],'k','FaceAlpha',.5); 
        hold off, axp = [min(0,min(min(-P))) max(0,max(max(-P)))];
        axis([ax axp]), title('level set function')
        drawnow
    end
end

%=======================================================================
function result= isInside(x,y,xc,yc,s)
    if abs((y-yc)) <= s/2 && abs((x-xc)) <= s/2
        result = true;
    else
        result = false;
    end
%--------------------------------------------------------------------------
% computing the gradient
function grad = gradi(P,h)
% computes normal 
Px = (P(3:end,:)-P(1:end-2,:))/(2*h); Px = Px([1 1:end end],:);
Py = (P(:,3:end)-P(:,1:end-2))/(2*h); Py = Py(:,[1 1:end end]);
grad = sqrt((Px.^2+Py.^2).^2);
% grad = min(max(F,-1/h),1/h);

%--------------------------------------------------------------------------
% computing the phix and phiy godunov
function dP = FabsgradP(P,h,F,c)
if nargin<4
    c = 0;
if nargin<3
    F = 1; 
end
end
DxP = diff(P)/h;   DxmP = DxP([1 1:end],:); DxpP = DxP([1:end end],:);
DyP = diff(P')'/h; DymP = DyP(:,[1 1:end]); DypP = DyP(:,[1:end end]);
Np = sqrt(max(DxmP,0).^2+min(DxpP,0).^2+max(DymP,0).^2+min(DypP,0).^2);
Nm = sqrt(min(DxmP,0).^2+max(DxpP,0).^2+min(DymP,0).^2+max(DypP,0).^2);
dP = max(F,0).*(Np-c)+min(F,0).*(Nm-c);
%--------------------------------------------------------------------------
function F = curvature(P,h)
% computes curvature by central differences
Pxx = diff(P([1 1:end end],:),2)/h^2;
Pyy = diff(P(:,[1 1:end end])',2)'/h^2;
Px = (P(3:end,:)-P(1:end-2,:))/(2*h); Px = Px([1 1:end end],:);
Py = (P(:,3:end)-P(:,1:end-2))/(2*h); Py = Py(:,[1 1:end end]);
Pxy = (Px(:,3:end)-Px(:,1:end-2))/(2*h); Pxy = Pxy(:,[1 1:end end]);
F = (Pxx.*Py.^2-2*Px.*Py.*Pxy+Pyy.*Px.^2)./(Px.^2+Py.^2).^1.5;
F = min(max(F,-1/h),1/h);