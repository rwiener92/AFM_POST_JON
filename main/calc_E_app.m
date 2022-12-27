function [E_app,regimeChange,varargout] = calc_E_app(D,F,R,th,b,mode,plotOpt)
% Calculates E_app using blunted cone model.
% Inputs:
% D: Depth vector Nx1 in nm
% F: Force vector Nx1 in N/m
% b, R, th : AFM tip parameters
% mode: 'pointwise' or 'Hertz' for pointwise or single value of E returned
% plotOpt: 1 or 0 to plot and not plot the linear plot for Hertz analysis
%
% Outputs:
% E_app : Vector of depthwise E_app values if mode is 'pointwise' or single
% value if mode is 'Hertz'
% varargout{1} : r^2 value for the linear fit if 'Hertz' mode is used

% Previous versions of code only have 5 inputs. This makes it backwards
% compatible.
if nargin == 5
    mode = 'pointwise';
elseif nargin == 6
    plotOpt = 0;
end

regimeChange = 0;
if strcmp(mode,'pointwise')
    E_app = zeros(length(F),1);
    for j = 1:length(D)
        % The solution for the modulus uses a spherical geometry until
        % the contact radius surpasses the cylindrical radius
        % and a blunt-cone geometry once this point is reached:
        Dj = D(j);
        if  Dj <= b^2/R %spherical
            E_app(j) = F(j) ./ (8/3*sqrt(Dj^3*R));
        else %blunted cone
            % For this geometry, must solve a nonlinear 1-d eqn for contact
            % area.
            [a,flag] = get_contact_radius_Jon(b,Dj,R,th);
            tm1 = a*Dj;
            tm2 = a^2 / (2*tan(th));
            tm2 = tm2 * (pi/2 - asin(b/a));
            tm3 = a^3/(3*R);
            tm4 = b/(2*tan(th));
            tm4 = tm4 + ((a^2-b^2)/(3*R)); % 3*R originally missing ()'s
            tm4 = tm4 * sqrt(a^2-b^2);
            E_app(j) = F(j) ./ (4*(tm1-tm2-tm3+tm4));
            if regimeChange == 0
                regimeChange = j;
            end
        end
    end
elseif strcmp(mode,'Hertz')
    %make sure F is row column
    [r c] = size(F);
    if c > r
       F = F'; 
    end
    x_fit = zeros(size(F));
    for j = 1:length(D)
        % The solution for the modulus uses a spherical geometry until
        % the contact radius surpasses the cylindrical radius
        % and a blunt-cone geometry once this point is reached:
        Dj = D(j);
        if  Dj <= b^2/R %spherical
            
            x_fit(j) = (8/3*sqrt(Dj^3*R));
        else %blunted cone
            % For this geometry, must solve a nonlinear 1-d eqn for contact
            % area.
            [a,flag] = get_contact_radius_Jon(b,Dj,R,th);
            tm1 = a*Dj;
            tm2 = a^2 / (2*tan(th));
            tm2 = tm2 * (pi/2 - asin(b/a));
            tm3 = a^3/(3*R);
            tm4 = b/(2*tan(th));
            tm4 = tm4 + ((a^2-b^2)/(3*R)); % 3*R originally missing ()'s
            tm4 = tm4 * sqrt(a^2-b^2);
            x_fit(j) =  4*(tm1-tm2-tm3+tm4);
            if regimeChange == 0
                regimeChange = j;
            end
        end
    end
    % now fit the linearized equations, the slope will be E.
    
    % This method includes fits y = mx + b. the b term has no physical
    % meaning and is ignored.
    %p = polyfit(x_fit,F,1);
    %E_app = p(1);
    
   % This method fits y = mx.
   % F = xE
   % x^TF = x^TxE
   % b = Ax
   E_app = (transpose(x_fit)*x_fit) \ (transpose(x_fit)*F);
   
   % Calculate R^2
   f = E_app*x_fit;
   ymean = mean(F);
   SSres = sum((F-f).^2);
   SStot = sum((F-ymean).^2);
   rsq = 1 - SSres/SStot;
   varargout{1} = rsq;

   
   
   
   
   if plotOpt == 1
      figure
      subplot(2,1,1)
      plot(D,F,'-*'); title('Force Curve'); xlabel('Depth (nm)');ylabel('Force (nN)');
      set(gca,'fontsize',14)
      subplot(2,1,2)
      plot(x_fit,F,'*'); title('Linearized Hertz'); xlabel('x') ;ylabel('Force (nN)');
      hold on; plot(x_fit,x_fit*E_app,'r-');
      legend('Raw','Linear Fit','Location','best')
      s = 'r^2 = %1.3f';
      s = sprintf(s,rsq);
      text(.1,.8,s,'Units','normalized','fontsize',16)
      set(gca,'fontsize',18)
      
      
   end
else
    error('Unknown mode.');
end

end

