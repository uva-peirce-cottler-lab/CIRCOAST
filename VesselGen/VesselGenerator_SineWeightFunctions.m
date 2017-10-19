function kernel_fcn = VesselGenerator_SineWeightFunctions(kernel, scale)

% keyboard
scale = 700;
kernel = zeros(3,3);
% y = a*sine(b*x + c)

c = rand(numel(kernel),3);
for n=1:numel(kernel);
    [r, c] = ind2sub(size(kernel),n);
    rc_str = sprintf('%.f,%.f', r, c);
    fcns{n}= [sprintf('%.4f',c(n,1)) ...
        ' * sin(' sprintf('%.4f', 1/(c(n,2)*scale)) ' * x(' rc_str ') + ' ...
        sprintf('%.4f', c(n,3)) '*pi)'];
end

keyboard
% n = 1;

% kernel_fxn = @(x) [fcns{1}(x(1,1)), fcns{2}(x(1,2))]