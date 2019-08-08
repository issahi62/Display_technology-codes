fh = figure;
screenfig = 5; 
for ii = 1:screenfig
    subplot(3,2,ii)
    P{ii} = get(gca,'pos');
end
clf
F = findobj('type','figure');
for ii = 2:6
    ax = findobj(F(ii),'type','axes');
    set(ax,'parent',fh,'pos',P{7-ii})
    close(F(ii))
end