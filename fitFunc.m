function c_disp = fitFunc(x,y,fitmode)
% curve fitting

if nargin<3,
    ord = [2 1]; % quad_origin
elseif ischar(fitmode),
    switch fitmode,
        case 'quadratic', % quadratic curve, (3 parameters)
            ord=[2 1 0];
        case 'quad_origin', % quadratic curve, but constrained to origin (2 parameters)
            X = [x(:).*abs(x(:)) x(:)];
            fprintf('fitting: p = c2*Q*|Q| + c1*Q\n'); ord=[2 1];
        case 'squared', % quadratic curve, pure y=c2*x^2 relationship (1 parameter)
            X = [x(:).*abs(x(:))];
            fprintf('fitting: p = c2*Q*|Q| + c1*Q\n'); ord=[2];
        case 'linear', % linear curve, (2 parameters)
            X = [x(:) ones(size(x(:)))];
            fprintf('fitting: p = c1*Q + c0\n'); ord=[1 0];
        case 'lin_origin', % linear curve, but constrained to origin (1 parameter)
            X = [x(:)];
            fprintf('fitting: p = c1*Q\n'); ord=[1];
        otherwise, % same as quad_origin
            X = [x(:).*abs(x(:)) x(:)];
            fprintf('fitting: p = c2*Q*|Q| + c1*Q\n'); ord=[2 1];
    end
else
    ord = fitmode;
end

fprintf('fitting: p = ');
if find(ord==2),
    fprintf('c2*Q*|Q| + ');
end
if find(ord==1),
    fprintf('c1*Q + ');
end
if find(ord==0),
    fprintf('c0 + ');
end
fprintf('\b\b\n');

X = [x(:).*abs(x(:)) x(:) ones(size(x(:)))];
X = X(:,3-ord);
c = robustfit(X,y,[],[],'off');
c_disp=[0 0 0];
c_disp(3-ord) = c;
