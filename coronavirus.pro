pro coronavirus

;; data here
;; https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series
  
str = importascii('~/personal/coronavirus.txt',/header)
g1=where(str.all gt 0)
g2=where(str.us gt 0)

setdisp
!p.font = 0
file='coronavirus'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=8.5
plot,[1],/nodata,xr=[1,75],yr=[1,2e4],xs=1,ys=1,/ylog,xtit='Date',ytit='New Confirmed Cases',$
     charsize=1.5,title='New Confirmed Coronavirus Cases',$
     xminor=4,xticks=3,xtickv=[1,20,40,60],xtickn=['Jan 21','Feb 9','Feb 29','Mar 20']
;     xminor=4,xticks=4,xtickv=[1,20,40,60,80],xtickn=['Jan 21','Feb 9','Feb 29','Mar 20','Apr 9']

oplot,str[g1].num,str[g1].all,ps=8,co=250
x=findgen(100)+30
oplot,x,0.3*exp(x/5),linestyle=2,co=250
; shifted by +25 days
;dgray = fsc_color('dark gray',1)
;oplot,str[g1].num+25,str[g1].all,ps=8,co=205
;x=findgen(100)+30
;oplot,x+25,0.3*exp(x/5),linestyle=2,co=205
; US
oplot,str[g2].num,str[g2].us,ps=8,co=70,sym=1.5
x=findgen(100)+14
oplot,x+28,0.3*exp(x/3.5),linestyle=2,co=80

gray = fsc_color('gray',2)
oplot,[0,56.5],[1000,1000],linestyle=1
oplot,[56.5,56.5],[5,1000],linestyle=1
oplot,[56.5,56.5],[1,2.1],linestyle=1

xyouts,56.5,4,'3 days to',align=0.5,co=0,charsize=1.1
xyouts,56.5,3,'1000 new US',align=0.5,co=0,charsize=1.1
xyouts,56.5,2.3,'cases a day',align=0.5,co=0,charsize=1.1

oplot,[0,64.5],[10000,10000],linestyle=1
oplot,[64.5,64.5],[1,30],linestyle=1
oplot,[64.5,64.5],[100,10000],linestyle=1

xyouts,66,70,'11 days to',align=0.5,co=0,charsize=1.1
xyouts,66,50,'10,000 new US',align=0.5,co=0,charsize=1.1
xyouts,66,37,'cases a day',align=0.5,co=0,charsize=1.1

;oplot,[0,65],[1000,1000],linestyle=1
;oplot,[65,65],[5,1000],linestyle=1
;oplot,[65,65],[1,2.1],linestyle=1
;
;xyouts,65,4,'18 days to',align=0.5,co=0,charsize=1.1
;xyouts,65,3,'1000 new US',align=0.5,co=0,charsize=1.1
;xyouts,65,2.3,'cases a day',align=0.5,co=0,charsize=1.1
;
;oplot,[0,77],[10000,10000],linestyle=1
;oplot,[77,77],[1,30],linestyle=1
;oplot,[77,77],[100,10000],linestyle=1
;xyouts,77,70,'30 days to',align=0.5,co=0,charsize=1.1
;xyouts,77,50,'10,000 new US',align=0.5,co=0,charsize=1.1
;xyouts,77,37,'cases a day',align=0.5,co=0,charsize=1.1


legend_old,['World minus China','US'],textcolor=[250,70],pos=[1,8000],charsize=1.7,box=0
;legend_old,['World-China','   shifted by +25 days','US'],textcolor=[250,205,70],pos=[1,8000],charsize=1.7,box=0

;bindata,str[g1].num,str[g1].all,xbin1,ybin1,binsize=5,/med
;oplot,xbin1,ybin1,ps=1,co=200

ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell


x=findgen(100)+30   
new = 0.3*exp(x/5)
tot = total(new,/cum)
ind = first_el(where(tot gt 300e6))
; x=96 but we need to shift it by +25, so day 121
; today is day 48, so another 73 days or ~May 22

stop

end
