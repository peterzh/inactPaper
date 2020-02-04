function cm = PurpleWhiteGreen(n, gamma)
% PurpleWhiteGreen(n, gamma)
% colormap going from purple to white to green
% n is number of distinct colors
% gamma is a gamma correction term (use to emphasize small values)

if nargin<1; n = 100; end
if nargin<2; gamma = 0.6; end

linespace(
cm = ([n*ones(1,n), n:-1:0 ; ...
      0:n, n-1:-1:0; ...
      0:n, ones(1,n)*n]' / n).^gamma;
% colormap(cm);
