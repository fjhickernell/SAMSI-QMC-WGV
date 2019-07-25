%Printing in eps problem MWE
figure
plot([0 1], [1 1], '-') %horizontal line
axis([0 1 -1.2 1.2]) %but want the axes as given
set(gca,'Visible','off') %turned off, but box to stay the same size
print -depsc HorizLine.eps %figure is compressed to thin box
print -dpng HorizLine.png %figure retains shape, but has extra margin compared to eps
