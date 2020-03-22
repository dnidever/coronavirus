pro coronavirus

;; Data from WHO coronavirus daily situation reports
;; https://www.who.int/emergencies/diseases/novel-coronavirus-2019/situation-reports
  
str = importascii('coronavirus.txt',/header)

setdisp
!p.font = 0
file='coronavirus'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=8.5
xr = [1,max(str.num)+15]
yr = [1,max(str.all)*2]
plot,[1],/nodata,xr=xr,yr=yr,xs=1,ys=1,/ylog,xtit='Date',ytit='New Confirmed Cases',$
     charsize=1.5,title='New Confirmed Coronavirus Cases',$
     xminor=4,xticks=3,xtickv=[1,20,40,60],xtickn=['Jan 21','Feb 9','Feb 29','Mar 20']

;; World minus China
gdw = where(str.all gt 0 and str.num ge 30,ngdw)
wcoef = robust_poly_fit(str[gdw].num,alog10(str[gdw].all),1)
g1 = where(str.all gt 0)
oplot,str[g1].num,str[g1].all,ps=8,co=250
x = findgen(100)+28
oplot,x,10^poly(x,wcoef),co=250,linestyle=2
;; doubling time
wdouble = alog10(2)/wcoef[1]
xyouts,5,600,'Doubling Times:',align=0,charsize=1.7,co=0
xyouts,5,350,stringize(wdouble,ndec=1)+' days',align=0,charsize=1.7,co=250

;; US
gdus = where(str.us gt 0 and str.num ge 44,ngdus)
uscoef = robust_poly_fit(str[gdus].num,alog10(str[gdus].us),1)
g2 = where(str.us gt 0)
oplot,str[g2].num,str[g2].us,ps=8,co=70,sym=1.5
x = findgen(100)+44
oplot,x,10^poly(x,uscoef),co=80,linestyle=2
;; doubling time
usdouble = alog10(2)/uscoef[1]
xyouts,5,200,stringize(usdouble,ndec=1)+' days',align=0,charsize=1.7,co=70


;; lines for 1,000 cases
gray = fsc_color('gray',2)
oplot,[0,56.5],[1000,1000],linestyle=1
oplot,[56.5,56.5],[1,1000],linestyle=1
;oplot,[56.5,56.5],[5,1000],linestyle=1
;oplot,[56.5,56.5],[1,2.1],linestyle=1

;xyouts,56.5,4,'1 day to',align=0.5,co=0,charsize=1.1
;xyouts,56.5,3,'1000 new US',align=0.5,co=0,charsize=1.1
;xyouts,56.5,2.3,'cases a day',align=0.5,co=0,charsize=1.1

;; lines for 10,000 cases a new
g = first_el(where(10^poly(x,uscoef) ge 1e4,ng))
x1e4 = x[g[0]]
oplot,[0,x1e4],[10000,10000],linestyle=1
oplot,[0,0]+x1e4,[1,30],linestyle=1
oplot,[0,0]+x1e4,[100,10000],linestyle=1

ndays1e4 = x1e4-max(str.num)
xyouts,64,70,strtrim(long(ndays1e4),2)+' days to',align=0.5,co=0,charsize=1.1
xyouts,64,50,'10,000 new US',align=0.5,co=0,charsize=1.1
xyouts,64,37,'cases a day',align=0.5,co=0,charsize=1.1

;; Today's date
jd = systime(/julian)
caldat,jd,month,day,year,hour,minute,second
months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
datestr = months[month-1]+' '+strtrim(long(day),2)+', '+strtrim(long(year),2)
xyouts,5,yr[1]*0.5,datestr,align=0,charsize=1.8,co=0

legend_old,['World minus China','US'],textcolor=[250,70],charsize=1.7,box=0,pos=[2,9000]
;legend_old,[datestr,'','World minus China','US'],textcolor=[0,255,250,70],charsize=1.7,box=0,/top,/left
;pos=[1,8000]

ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell

stop

end
