;*************************************************
; CT-NRTvs.aircraft.usepbl.ncl
; based on overlay_1.ncl
;
; Concepts illustrated:
;   - Overlaying line contours on filled contours
;   - Explicitly setting contour levels
;   - Selecting a different color map for each contour plot
;
;*************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

;;; Read WW3 output depth netcdf file, in order to generate ASCII text output
;;; needed for running scrip_coawst and generating script netcdf mapping file.


;  WW3file="../../../../../wave/2014/bulk/ww3.bulk.vi_hr.201407_dpt.nc"
;;;  WW3file="./ww3.bulk.at_2km_new.201709_dpt.nc"
  WW3file="./ww3.bulk.at_2km_large.201709_dpt.nc"
  a = addfile(WW3file,"r")
  xcoords=a->longitude
  ycoords=a->latitude
  bath=a->dpt
  delete(a)

  ww3_dimarray = dimsizes(xcoords)
  nx_ww3 = ww3_dimarray(0)
  delete(ww3_dimarray)

  ww3_dimarray = dimsizes(ycoords)
  ny_ww3 = ww3_dimarray(0)
  delete(ww3_dimarray)

;; iprint = 1:  print x coordinate file
;; iprint = 2:  print y coordinate file
;; iprint = 3:  print bathymetry file

iprint = 3 

  do iy=0,ny_ww3-1
     do ix=0,nx_ww3-1
       if( iprint .eq. 1 ) then 
          print(xcoords(ix)-360.0)
       else if (iprint .eq. 2 ) then
         print(ycoords(iy)+" ")
       else if (iprint .eq. 3 ) then

         if(ismissing(bath(0,iy,ix))) then
              print("9999")
         else if(bath(0,iy,ix).gt.1.0e20) then
             print("9999")
         else
             print(bath(0,iy,ix)+" ")
         end if
         end if
       end if
       end if
       end if
     end do
  end do



end
   
  
