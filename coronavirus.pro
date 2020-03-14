pro coronavirus

;; Data from WHO coronavirus daily situation reports
;; https://www.who.int/emergencies/diseases/novel-coronavirus-2019/situation-reports
  
str = importascii('coronavirus.txt',/header)

setdisp
!p.font = 0
file='coronavirus'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=8.5
plot,[1],/nodata,xr=[1,75],yr=[1,2e4],xs=1,ys=1,/ylog,xtit='Date',ytit='New Confirmed Cases',$
     charsize=1.5,title='New Confirmed Coronavirus Cases',$
     xminor=4,xticks=3,xtickv=[1,20,40,60],xtickn=['Jan 21','Feb 9','Feb 29','Mar 20']

;; World minus China
gdw = where(str.all gt 0 and str.num ge 30,ngdw)
wcoef = robust_poly_fit(str[gdw].num,alog10(str[gdw].all),1)
g1 = where(str.all gt 0)
oplot,str[g1].num,str[g1].all,ps=8,co=250
x = findgen(100)+28
oplot,x,10^poly(x,wcoef),co=250,linestyle=2

;; US
gdus = where(str.us gt 0 and str.num ge 44,ngdus)
uscoef = robust_poly_fit(str[gdus].num,alog10(str[gdus].us),1)
g2 = where(str.us gt 0)
oplot,str[g2].num,str[g2].us,ps=8,co=70,sym=1.5
x = findgen(100)+44
oplot,x,10^poly(x,uscoef),co=80,linestyle=2

;; lines for 1,000 cases
gray = fsc_color('gray',2)
oplot,[0,56.5],[1000,1000],linestyle=1
oplot,[56.5,56.5],[5,1000],linestyle=1
oplot,[56.5,56.5],[1,2.1],linestyle=1

xyouts,56.5,4,'3 days to',align=0.5,co=0,charsize=1.1
xyouts,56.5,3,'1000 new US',align=0.5,co=0,charsize=1.1
xyouts,56.5,2.3,'cases a day',align=0.5,co=0,charsize=1.1

;; lines for 10,000 cases a new
oplot,[0,64.5],[10000,10000],linestyle=1
oplot,[64.5,64.5],[1,30],linestyle=1
oplot,[64.5,64.5],[100,10000],linestyle=1

xyouts,66,70,'11 days to',align=0.5,co=0,charsize=1.1
xyouts,66,50,'10,000 new US',align=0.5,co=0,charsize=1.1
xyouts,66,37,'cases a day',align=0.5,co=0,charsize=1.1

legend_old,['World minus China','US'],textcolor=[250,70],pos=[1,8000],charsize=1.7,box=0

ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell

stop

end
