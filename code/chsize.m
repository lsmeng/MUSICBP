function chsize(N)

set(gca,'FontSize',N)
h_xlabel = get(gca,'XLabel');
set(h_xlabel,'FontSize',N); 
h_ylabel = get(gca,'YLabel');
set(h_ylabel,'FontSize',N); 
h_title = get(gca,'Title');
set(h_title,'FontSize',N); 