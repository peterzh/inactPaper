function fig2Pdf(f,filepath)

set(f,'Units','Inches','renderer','painters');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(f,filepath,'-dpdf','-r0')


end