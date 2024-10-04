function [g, badg] = numgrad(fcn,x,varargin)
global delta_
% function [g badg] = numgrad(fcn,xvarargin)
%
% delta = delta_;
delta = 1e-7;
delta = 1e-13;
% delta = 1e-6;
% delta=1e-2;

n=length(x);

tvec=delta*eye(n);
% tvec = delta*(max(abs(x),1e-8*ones(n,1)));

% g=zeros(n,1);
%--------------------old way to deal with variable # of P's--------------
%tailstr = ')';
%stailstr = [];
%for i=nargin-2:-1:1
%   tailstr=[ ',P' num2str(i)  tailstr];
%   stailstr=[' P' num2str(i) stailstr];
%end
%f0 = eval([fcn '(x' tailstr]); % Is there a way not to do this?
%---------------------------------------------------------------^yes
f0 = eval([fcn '(x,varargin{:})']);
g=zeros(length(f0),n);

% disp(' first fcn in numgrad.m ------------------')
%home
% disp('numgrad.m is working. ----') % Jiinil on 9/5/95
% sizex=size(x),sizetvec=size(tvec),x,    % Jinill on 9/6/95
badg=0;
for i=1:n
    scale=1; % originally 1
    % i,tveci=tvec(:,i)% ,plus=x+scale*tvec(:,i) % Jinill Kim on 9/6/95
    if size(x,1)>size(x,2)
        tvecv=tvec(i,:);
    else
        tvecv=tvec(:,i);
    end
    g0 = (eval([fcn '(x+scale*tvecv'', varargin{:})']) - f0) ...
       /(scale*delta);
%     f1  = (eval([fcn '(x+scale*tvecv'', varargin{:})']));
%     f2  = (eval([fcn '(x-scale*tvecv'', varargin{:})']));
%     g0  = (f1-f2)/(2*scale*delta);
    
    
    % disp(' fcn in the i=1:n loop of numgrad.m ------------------')% Jinill 9/6/95
    % disp('          and i is')               % Jinill
    % i                         % Jinill
    % fprintf('Gradient w.r.t. %3d: %10g\n',i,g0) %see below Jinill 9/6/95
    % -------------------------- special code to essentially quit here
    % absg0=abs(g0) % Jinill on 9/6/95
    
    % original code
    %    if abs(g0)< 1e15
    %       g(i)=g0;
    %       % disp('good gradient') % Jinill Kim
    %    else
    %       disp('bad gradient ------------------------') % Jinill Kim
    %       % fprintf('Gradient w.r.t. %3d: %10g\n',i,g0) %see above
    %       g(i)=0;
    %       badg=1;
    %       % return
    %       % can return here to save time if the gradient will never be
    %       % used when badg returns as true.
    %    end
    
    % from FS

    g(:,i) = g0;
    if sqrt(norm(g0)) > 1e20
        
        gb = (f0 - eval([fcn '(x-scale*tvecv'', varargin{:})']))/(scale*delta);
        if (gb > 0) && (abs(gb) < 1e15)
            
            g(i) = gb;
        else
            g(i) = 0;
            badg = 1;
        end
    end
end


end
%-------------------------------------------------------------
%     if g0 > 0
%        sided=2;
%        g1 = -(eval([fcn '(x-scale*tvec(:,i)''' tailstr]) - f0) ...
%           /(scale*delta);
%        if g1<0
%           scale = scale/10;
%        else
%           break
%        end
%     else
%        sided=1;
%        break
%     end
%  end
%  if sided==1
%     g(i)=g0;
%  else
%     if (g0<1e20)
%        if (g1>-1e20)
%           g(i)=(g0+g1)/2;
%        else
%           g(i)=0;
%           badg=1;
%           disp( ['Banging against wall, parameter ' int2str(i)] );
%        end
%     else
%        if g1>-1e20
%           if g1<0
%              g(i)=0;
%              badg=1;
%              disp( ['Banging against wall, parameter ' int2str(i)] );
%           else
%              g(i)=g1;
%           end
%        else
%           g(i)=0;
%           badg=1;
%           disp(['Valley around parameter ' int2str(i)])
%        end
%     end
%  end
%end
%save g.dat g x f0
%eval(['save g g x f0 ' stailstr]);
