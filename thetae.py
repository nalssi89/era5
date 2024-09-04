import matplotlib
matplotlib.use('Agg')  # GUI 없이 이미지 생성

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# 상당온위 계산 함수 (간단한 근사식 사용)
def calc_equivalent_potential_temperature(temperature, relative_humidity, pressure):
    # 온도를 켈빈으로 변환
    T = temperature + 273.15
    # 포화수증기압 계산 (Bolton 1980)
    es = 6.112 * np.exp(17.67 * (T - 273.15) / (T - 29.65))
    # 혼합비 계산
    r = 0.622 * (relative_humidity / 100) * es / (pressure - es)
    # 상당온위 계산 (Bolton 1980의 근사식)
    theta_e = T * (1000 / pressure) ** 0.2854 * (1 + 0.28 * r) * np.exp((3376 / T - 2.54) * r * (1 + 0.81 * r))
    return theta_e

# netCDF 파일 읽기
ds = xr.open_dataset('era5_prs.nc')

# 필요한 변수 추출
temp = ds['t'].squeeze().values - 273.15  # 온도 (°C)
rh = ds['r'].squeeze().values    # 상대습도 (%)
pressure = 850  # 압력 (hPa)

# 상당온위 계산
equiv_potential_temp = calc_equivalent_potential_temperature(temp, rh, pressure)

# 그림 설정
plt.figure(figsize=(15, 10))

# 지도 투영법 설정 (여기서는 PlateCarree 투영법 사용)
proj = ccrs.PlateCarree()

# 지도 생성
ax = plt.axes(projection=proj)

# 해안선 추가
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS)

# 상당온위 등치선 그리기 (3K 간격)
min_level = int(np.floor(equiv_potential_temp.min()))
max_level = int(np.ceil(equiv_potential_temp.max()))
levels = np.arange(min_level, max_level, 3)

# 무지개 색상으로 등치선 그리기
contour = ax.contourf(ds.longitude, ds.latitude, equiv_potential_temp, 
                      levels=levels, cmap='rainbow', transform=proj, extend='both')

# 등치선 레이블 추가
contour_lines = ax.contour(ds.longitude, ds.latitude, equiv_potential_temp, 
                           levels=levels, colors='k', linewidths=0.5, transform=proj)
ax.clabel(contour_lines, inline=True, fontsize=8, fmt='%d')

# 그림 꾸미기
plt.title('Equivalent Potential Temperature at 850 hPa')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# 경위도 그리드 추가
ax.gridlines(draw_labels=True, linestyle='--', alpha=0.6)

# 컬러바 추가
cbar = plt.colorbar(contour, ax=ax, orientation='horizontal', pad=0.05)
cbar.set_label('Equivalent Potential Temperature (K)')

# 지도 범위 설정 (동아시아 지역으로 설정)
ax.set_extent([100, 150, 10, 50], crs=proj)

# 그림 저장
plt.savefig('equiv_potential_temp_850hPa_map.png', dpi=300, bbox_inches='tight')

print("그림이 성공적으로 저장되었습니다: equiv_potential_temp_850hPa_map.png")