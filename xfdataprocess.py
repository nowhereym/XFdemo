import os
import numpy as np
import netCDF4 as nc
import logging
import matplotlib.pyplot as plt
from pykrige import OrdinaryKriging
from scipy.interpolate import griddata, NearestNDInterpolator
from scipy.spatial import cKDTree
from sklearn.metrics import mean_squared_error, mean_absolute_error,r2_score
from scipy.spatial import KDTree
from scipy.stats import pearsonr


# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def bilinear_interpolation( data_9km, target_lon, target_lat, variable):
    """双线性插值"""
    points = np.column_stack([data_9km['lon'].flatten(), data_9km['lat'].flatten()])
    values = data_9km[variable].flatten()

    # 移除NaN值
    mask = ~np.isnan(values)
    points = points[mask]
    values = values[mask]

    # 使用scipy的griddata进行双线性插值
    interpolated = griddata(points, values, (target_lon, target_lat),
                            method='linear', fill_value=np.nan)
    return interpolated

def load_netcdf_data(file_path):
    """加载NetCDF格式的气象数据"""
    logger.info(f"正在加载NetCDF文件: {file_path}")
    ds = nc.Dataset(file_path)
    # 打印数据集信息用于调试
    logger.info(f"数据集变量: {list(ds.variables.keys())}")
    return ds

def tiqvdata(ds):
    """从数据集中提取气象变量，处理[H, W, T]维度顺序"""
    lon=ds['longitude']
    lat=ds['latitude']
    # 使用用户指定的变量名
    variable_mapping = {
        't2m': ['t2m'],
        'v10': ['v10'],
        'u10': ['u10'],
        'tp6h': ['tp6h'],
        'ssr6h': ['ssr6h']
    }
    extracted_data = {}
    for target_var, possible_names in variable_mapping.items():
        for name in possible_names:
            if name in ds.variables:
                data = ds[name]
                logger.info(f"找到变量 {name}, 形状: {data.shape}")
                # 处理不同的维度顺序
                data = np.transpose(data, (2, 0, 1))
                extracted_data[target_var] = data
                logger.info(f"提取变量 {target_var}, 最终形状: {data.shape}")
        # 记录所有提取的变量
    logger.info(f"成功提取的变量: {list(extracted_data.keys())}")
    return extracted_data,lat,lon

def idw_interpolation(tree, values, target_points, k=4, power=2):
    """
    反距离加权插值

    参数:
    - tree: 源数据点的cKDTree
    - values: 源数据点的值
    - target_points: 目标点坐标
    - k: 使用的最近邻点数
    - power: 距离的幂次

    返回:
    - 插值结果
    """
    # 查询k个最近邻的距离和索引
    distances, indices = tree.query(target_points, k=k)

    # 避免除零错误，将零距离替换为一个小值
    distances = np.where(distances == 0, 1e-12, distances)

    # 计算权重 (1 / distance^power)
    weights = 1.0 / (distances ** power)

    # 对每个目标点，计算加权平均值
    weighted_sum = np.sum(weights * values[indices], axis=1)
    sum_weights = np.sum(weights, axis=1)

    # 避免除零错误
    sum_weights = np.where(sum_weights == 0, 1e-12, sum_weights)

    return weighted_sum / sum_weights

def ordinary_kriging_interpolation(lon_src, lat_src, values_src, lon_target, lat_target,variogram_model='linear', nlags=6, weight=True):
    """
    普通克里金插值
    """
    # 移除NaN值
    mask = ~np.isnan(values_src)
    lon_src = lon_src[mask]
    lat_src = lat_src[mask]
    values_src = values_src[mask]

    # 创建普通克里金对象
    OK = OrdinaryKriging(
        lon_src,
        lat_src,
        values_src,
        variogram_model=variogram_model,
        nlags=nlags,
        weight=weight,
        verbose=False,
        enable_plotting=False
    )

    # 执行插值
    z, ss = OK.execute('grid', lon_target, lat_target)
    return z

def analyse_extract_meteorological_data(extracted_data,lat,lon,var,era5,elon,elat):
    """
     三种插值：将9km分辨率数据插值到1km分辨率
     参数:
     - data_9km: 形状为 (265, 423, 723) 的气象数据
     返回:
     - data_1km: 形状为 (265, 3807, 6507) 的1km分辨率数据
     误差分布空间图
     """
    # 获取原始数据维度
    # 新疆:lon[11:268],lat[189:357]
    ds=extracted_data[var][:,189:357,11:268]
    n_time, n_lat, n_lon = ds.shape
    print(n_time, n_lat, n_lon)
    # 计算1km分辨率的维度 (9倍于9km分辨率)
    n_lat_1km = n_lat * 9
    n_lon_1km = n_lon * 9
    n_lat_25km = int(n_lat_1km/25)
    n_lon_25km = int(n_lon_1km/25)
    print(f"原始数据形状: {ds.shape}")
    print(f"目标数据形状: ({n_time}, {n_lat_1km}, {n_lon_1km})")
    print(f"25km目标数据形状: ({n_time}, {n_lat_25km}, {n_lon_25km})")

    # 确保经纬度是1D数组
    if lon.ndim == 2:
        lon_1d = lon[0, :]  # 取第一行
    else:
        lon_1d = lon

    if lat.ndim == 2:
        lat_1d = lat[:, 0]  # 取第一列
    else:
        lat_1d = lat
    #era5经纬度

    lon_e=elon[:]
    lat_e=elat[:]
    # 创建1km目标网格
    lon_1km = np.linspace(np.min(lon_1d), np.max(lon_1d), n_lon_1km)
    lat_1km = np.linspace(np.min(lat_1d), np.max(lat_1d), n_lat_1km)
    lon_1km_grid, lat_1km_grid = np.meshgrid(lon_1km, lat_1km)
    # 创建25km目标网格
    lon_25km = np.linspace(np.min(lon_e), np.max(lon_e), n_lon_25km)
    lat_25km = np.linspace(np.min(lat_e), np.max(lat_e), n_lat_25km)
    lon_25km_grid, lat_25km_grid = np.meshgrid(lon_25km, lat_25km)
    # 创建原始网格点
    lon_9km_grid, lat_9km_grid = np.meshgrid(lon_1d, lat_1d)
    points = np.column_stack([lon_9km_grid.ravel(), lat_9km_grid.ravel()])
    target_points = np.column_stack([lon_1km_grid.ravel(), lat_1km_grid.ravel()])
    # 构建原始数据点的KD树用于快速最近邻搜索
    source_points = np.column_stack([lon_9km_grid.ravel(), lat_9km_grid.ravel()])
    tree = cKDTree(source_points)
    # 初始化结果数组
    #data_1km = np.zeros((n_time, n_lat_1km, n_lon_1km))

    # 对每个时间步进行插值
    for t in range(1):
        print(f"处理时间步 {t + 1}/{n_time}")
        # 获取当前时间步的数据
        values = ds[t, :, :].ravel()
        # 执行双线性插值
        interpolated_linear = griddata(
            points,
            values,
            (lon_1km_grid, lat_1km_grid),
            method='linear',
            fill_value=np.nan
        )
        # 将1km数据降尺度为25km
        print("开始降尺度到25km...")
        downsampled_25km = downsample_to_25km(
            interpolated_linear,
            lon_1km_grid,
            lat_1km_grid,
            lon_25km_grid,
            lat_25km_grid,
            method='mean'
        )
        print(downsampled_25km.shape,era5.shape)
        #era5plot_simple_comparison(np.flipud(downsampled_25km),era5[t],np.flipud(downsampled_25km)-era5[t,:60,:92],var_name=var)
        # 计算基本评估指标
        rmse, mae, correlation = calculate_basic_metrics(
            era5[t,:60,:92],
            downsampled_25km
        )
        # 打印结果
        print_basic_metrics(rmse, mae, correlation,var)

        # 使用反距离加权插值
        interpolated_idw = idw_interpolation(tree, values, target_points, k=4, power=2).reshape(lon_1km_grid.shape)
        # 内存不够
        # interpolated_kriging= ordinary_kriging_interpolation(
        #     lon_9km_grid.ravel(), lat_9km_grid.ravel(), values,
        #     lon_1km, lat_1km)
        interpolator = NearestNDInterpolator(list(zip(lon_9km_grid.ravel(), lat_9km_grid.ravel())),values)
        interpolated_2d = interpolator(lon_1km_grid, lat_1km_grid)
        # 执行降采样评估
        downsampled_9km, errors, error_stats = downsample_evaluation_simple(
            ds[t, :, :], interpolated_2d,
            lat_9km_grid, lon_9km_grid,  # 9km网格
            lat_1km_grid, lon_1km_grid,  # 1km网格
            var)
        ##绘图
        #plot_three_images(ds[t, :, :],interpolated_linear,interpolated_idw,interpolated_2d)
        # 随机点验证
        a,b=random_point_validation(ds[t, :, :],interpolated_linear,lon_9km_grid,lat_9km_grid,lon_1km_grid,lat_1km_grid)
        print(a,b)
        #np.save(f'data_1km/{var}-{t}',interpolated)

#反插值方法找对应点
# def random_point_validation(data_9km, data_1km, lon_9km, lat_9km, lon_1km, lat_1km, n_samples=1):
#     """
#     随机点验证（自一致性检验）
#
#     参数:
#     - data_9km: 原始9km分辨率数据 (时间, 纬度, 经度)
#     - data_1km: 插值后的1km分辨率数据 (时间, 纬度, 经度)
#     - lon_9km, lat_9km: 原始9km网格坐标
#     - lon_1km, lat_1km: 插值后的1km网格坐标
#     - n_samples: 随机抽样点数
#
#     返回:
#     - metrics: 包含各种评估指标的字典
#     - validation_points: 验证点的详细信息
#     """
#     # 选择第一个时间步进行验证
#     #time_idx = 0
#     data_9km_single = data_9km[ :, :]
#     data_1km_single = data_1km[ :, :]
#
#     # 确保经纬度是2D网格
#     if lon_9km.ndim == 1 and lat_9km.ndim == 1:
#         lon_9km_grid, lat_9km_grid = np.meshgrid(lon_9km, lat_9km)
#     else:
#         lon_9km_grid, lat_9km_grid = lon_9km, lat_9km
#
#     if lon_1km.ndim == 1 and lat_1km.ndim == 1:
#         lon_1km_grid, lat_1km_grid = np.meshgrid(lon_1km, lat_1km)
#     else:
#         lon_1km_grid, lat_1km_grid = lon_1km, lat_1km
#
#     # 随机选择验证点
#     n_lat, n_lon = data_9km_single.shape
#     total_points = n_lat * n_lon
#
#     # 确保抽样点数不超过总点数
#     n_samples = min(n_samples, total_points)
#
#     # 随机选择点索引
#     random_indices = np.random.choice(total_points, n_samples, replace=False)
#     lat_indices, lon_indices = np.unravel_index(random_indices, (n_lat, n_lon))
#
#     # 获取验证点的原始值和坐标
#     original_values = data_9km_single[lat_indices, lon_indices]
#     validation_lons = lon_9km_grid[lat_indices, lon_indices]
#     validation_lats = lat_9km_grid[lat_indices, lon_indices]
#
#     # 从1km插值结果中获取验证点的值
#     # 使用双线性插值从1km网格获取验证点的值
#     points_1km = np.column_stack([lon_1km_grid.ravel(), lat_1km_grid.ravel()])
#     values_1km = data_1km_single.ravel()
#
#     interpolated_values = griddata(
#         points_1km, values_1km,
#         (validation_lons, validation_lats),
#         method='linear',
#         fill_value=np.nan
#     )
#     # 移除NaN值
#     original_values_valid = original_values
#     interpolated_values_valid = interpolated_values
#     validation_lons_valid = validation_lons
#     validation_lats_valid = validation_lats
#
#     # 计算评估指标
#     rmse = np.sqrt(mean_squared_error(original_values_valid, interpolated_values_valid))
#     mae = mean_absolute_error(original_values_valid, interpolated_values_valid)
#     r2 = r2_score(original_values_valid, interpolated_values_valid)
#
#     # 计算偏差和标准差
#     bias = np.mean(interpolated_values_valid - original_values_valid)
#     std_dev = np.std(interpolated_values_valid - original_values_valid)
#
#     # 计算平均绝对百分比误差 (MAPE)
#     mape = np.mean(np.abs((interpolated_values_valid - original_values_valid) / original_values_valid)) * 100
#
#     metrics = {
#         'rmse': rmse,
#         'mae': mae,
#         'r2': r2,
#         'bias': bias,
#         'std_dev': std_dev,
#         'mape': mape,
#         'n_valid_points': len(original_values_valid)
#     }
#
#     validation_points = {
#         'lons': validation_lons_valid,
#         'lats': validation_lats_valid,
#         'original_values': original_values_valid,
#         'interpolated_values': interpolated_values_valid
#     }
#
#     return metrics, validation_points

def random_point_validation(data_9km, data_1km, lon_9km, lat_9km, lon_1km, lat_1km, n_samples=3):
    """
    随机点验证（自一致性检验）

    参数:
    - data_9km: 原始9km分辨率数据 (时间, 纬度, 经度)
    - data_1km: 插值后的1km分辨率数据 (时间, 纬度, 经度)
    - lon_9km, lat_9km: 原始9km网格坐标
    - lon_1km, lat_1km: 插值后的1km网格坐标
    - n_samples: 随机抽样点数

    返回:
    - metrics: 包含各种评估指标的字典
    - validation_points: 验证点的详细信息
    """
    # 选择第一个时间步进行验证
    # time_idx = 0
    data_9km_single = data_9km[:, :]
    data_1km_single = data_1km[:, :]

    # 确保经纬度是2D网格
    if lon_9km.ndim == 1 and lat_9km.ndim == 1:
        lon_9km_grid, lat_9km_grid = np.meshgrid(lon_9km, lat_9km)
    else:
        lon_9km_grid, lat_9km_grid = lon_9km, lat_9km

    if lon_1km.ndim == 1 and lat_1km.ndim == 1:
        lon_1km_grid, lat_1km_grid = np.meshgrid(lon_1km, lat_1km)
    else:
        lon_1km_grid, lat_1km_grid = lon_1km, lat_1km

    # 随机选择验证点
    n_lat, n_lon = data_9km_single.shape
    total_points = n_lat * n_lon

    # 确保抽样点数不超过总点数
    n_samples = min(n_samples, total_points)

    # 随机选择点索引
    random_indices = np.random.choice(total_points, n_samples, replace=False)
    lat_indices, lon_indices = np.unravel_index(random_indices, (n_lat, n_lon))

    # 获取验证点的原始值和坐标
    original_values = data_9km_single[lat_indices, lon_indices]
    validation_lons = lon_9km_grid[lat_indices, lon_indices]
    validation_lats = lat_9km_grid[lat_indices, lon_indices]
    print(original_values, validation_lats, validation_lons)

    # 使用KDTree最近邻方法从1km网格获取验证点的值
    # 准备1km网格点坐标
    points_1km = np.column_stack([lon_1km_grid.ravel(), lat_1km_grid.ravel()])
    values_1km = data_1km_single.ravel()

    # 创建KD树
    tree = KDTree(points_1km)

    # 查询每个验证点的最近1km网格点
    validation_points_xy = np.column_stack([validation_lons, validation_lats])
    distances, indices = tree.query(validation_points_xy)
    interpolated_values = values_1km[indices]

    print(f"最大匹配距离: {np.max(distances):.6f} 度")
    print(interpolated_values.shape)

    # 移除NaN值
    valid_mask = ~np.isnan(interpolated_values)
    original_values_valid = original_values[valid_mask]
    interpolated_values_valid = interpolated_values[valid_mask]
    validation_lons_valid = validation_lons[valid_mask]
    validation_lats_valid = validation_lats[valid_mask]

    # 计算评估指标
    rmse = np.sqrt(mean_squared_error(original_values_valid, interpolated_values_valid))
    mae = mean_absolute_error(original_values_valid, interpolated_values_valid)
    r2 = r2_score(original_values_valid, interpolated_values_valid)

    # 计算偏差和标准差
    bias = np.mean(interpolated_values_valid - original_values_valid)
    std_dev = np.std(interpolated_values_valid - original_values_valid)

    # 计算平均绝对百分比误差 (MAPE)
    mape = np.mean(np.abs((interpolated_values_valid - original_values_valid) / original_values_valid)) * 100

    metrics = {
        'rmse': rmse,
        'mae': mae,
        'r2': r2,
        'bias': bias,
        'std_dev': std_dev,
        'mape': mape,
        'n_valid_points': len(original_values_valid),
        'max_distance': np.max(distances)  # 添加最大匹配距离作为参考
    }

    validation_points = {
        'lons': validation_lons_valid,
        'lats': validation_lats_valid,
        'original_values': original_values_valid,
        'interpolated_values': interpolated_values_valid,
        'distances': distances[valid_mask]  # 添加匹配距离信息
    }

    return metrics, validation_points


def downsample_to_25km(interpolated_1km, lon_1km_grid, lat_1km_grid, lon_25km_grid, lat_25km_grid, method='mean'):
    """
    将1km分辨率数据降尺度到25km分辨率

    参数:
    - interpolated_1km: 1km分辨率数据 (2D数组)
    - lon_1km_grid: 1km网格经度坐标 (2D数组)
    - lat_1km_grid: 1km网格纬度坐标 (2D数组)
    - lon_25km_grid: 25km网格经度坐标 (2D数组)
    - lat_25km_grid: 25km网格纬度坐标 (2D数组)
    - method: 降尺度方法，可选 'mean', 'median', 'max', 'min'

    返回:
    - downsampled_25km: 25km分辨率数据 (2D数组)
    """

    # 获取25km网格的形状
    n_lat_25km, n_lon_25km = lon_25km_grid.shape

    # 初始化25km结果数组
    downsampled_25km = np.full((n_lat_25km, n_lon_25km), np.nan)

    # 计算1km网格的经纬度间隔
    d_lat_1km = np.abs(lat_1km_grid[1, 0] - lat_1km_grid[0, 0])
    d_lon_1km = np.abs(lon_1km_grid[0, 1] - lon_1km_grid[0, 0])

    # 计算25km网格的经纬度间隔
    d_lat_25km = np.abs(lat_25km_grid[1, 0] - lat_25km_grid[0, 0])
    d_lon_25km = np.abs(lon_25km_grid[0, 1] - lon_25km_grid[0, 0])

    # 计算每个25km网格包含的1km网格数
    lat_ratio = int(d_lat_25km / d_lat_1km)
    lon_ratio = int(d_lon_25km / d_lon_1km)

    print(f"降尺度比例: 纬度方向 {lat_ratio}:1, 经度方向 {lon_ratio}:1")

    # 对每个25km网格进行降尺度
    for i in range(n_lat_25km):
        for j in range(n_lon_25km):
            # 计算当前25km网格的中心坐标
            center_lat = lat_25km_grid[i, j]
            center_lon = lon_25km_grid[i, j]

            # 计算在1km网格中对应的范围
            lat_min = center_lat - d_lat_25km / 2
            lat_max = center_lat + d_lat_25km / 2
            lon_min = center_lon - d_lon_25km / 2
            lon_max = center_lon + d_lon_25km / 2

            # 找到在范围内的1km网格点
            mask = ((lat_1km_grid >= lat_min) & (lat_1km_grid <= lat_max) &
                    (lon_1km_grid >= lon_min) & (lon_1km_grid <= lon_max))

            # 提取范围内的1km数据
            values_in_range = interpolated_1km[mask]

            # 根据选择的方法进行聚合
            if len(values_in_range) > 0 and np.any(~np.isnan(values_in_range)):
                if method == 'mean':
                    downsampled_25km[i, j] = np.nanmean(values_in_range)
                elif method == 'median':
                    downsampled_25km[i, j] = np.nanmedian(values_in_range)
                elif method == 'max':
                    downsampled_25km[i, j] = np.nanmax(values_in_range)
                elif method == 'min':
                    downsampled_25km[i, j] = np.nanmin(values_in_range)
                else:
                    downsampled_25km[i, j] = np.nanmean(values_in_range)

    return downsampled_25km

def downsample_evaluation_simple(original_9km, interpolated_1km, lat_9km, lon_9km,
                                 lat_1km, lon_1km, var_name):
    """
    简化的降采样评估
    """
    # 获取原始9km数据的维度
    n_lat_9km, n_lon_9km = original_9km.shape

    # 初始化降采样后的9km数据数组
    downsampled_9km = np.full((n_lat_9km, n_lon_9km), np.nan)

    # 计算降采样比例
    lat_ratio = int(interpolated_1km.shape[0] / n_lat_9km)
    lon_ratio = int(interpolated_1km.shape[1] / n_lon_9km)

    print(f"降采样比例: 纬度方向 {lat_ratio}:1, 经度方向 {lon_ratio}:1")

    # 执行降采样
    for i in range(n_lat_9km):
        for j in range(n_lon_9km):
            lat_start = i * lat_ratio
            lat_end = min((i + 1) * lat_ratio, interpolated_1km.shape[0])
            lon_start = j * lon_ratio
            lon_end = min((j + 1) * lon_ratio, interpolated_1km.shape[1])

            block = interpolated_1km[lat_start:lat_end, lon_start:lon_end]

            if np.any(~np.isnan(block)):
                downsampled_9km[i, j] = np.nanmean(block)

    # 计算误差
    errors = downsampled_9km - original_9km

    # 计算误差统计
    valid_errors = errors[~np.isnan(errors)]
    error_stats = {
        'mean_error': np.mean(valid_errors),
        'std_error': np.std(valid_errors),
        'rmse': np.sqrt(np.mean(valid_errors ** 2)),
        'mae': np.mean(np.abs(valid_errors)),
        'max_positive': np.max(valid_errors),
        'max_negative': np.min(valid_errors),
        'correlation': np.corrcoef(original_9km[~np.isnan(original_9km)],
                                   downsampled_9km[~np.isnan(downsampled_9km)])[0, 1]
    }

    # 打印误差统计
    print(f"\n=== {var_name}误差统计 ===")
    for key, value in error_stats.items():
        print(f"{key}: {value:.4f}")

    # 绘制简化图表
    plot_simple_comparison(original_9km, downsampled_9km, errors, var_name)

    return downsampled_9km, errors, error_stats

def era5plot_simple_comparison(original_25km, era5_25km, errors, var_name):
    """
    绘制简化比较图表
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 10))

    # 原始数据
    im1 = axes[0].imshow(original_25km, cmap='viridis', origin='lower')
    axes[0].set_title(f'{var_name} - 1km-to-25km')
    #plt.colorbar(im1, ax=axes[0])

    # 降采样数据
    im2 = axes[1].imshow(era5_25km, cmap='viridis', origin='lower')
    axes[1].set_title(f'{var_name} - era5-25km')
    #plt.colorbar(im2, ax=axes[1])

    # 误差分布
    im3 = axes[2].imshow(errors, cmap='RdBu_r', origin='lower')
    axes[2].set_title(f'{var_name} - difference')
    #plt.colorbar(im3, ax=axes[2])
    plt.tight_layout()
    plt.show()
    plt.savefig(f'era5dif-{var_name}')


def plot_simple_comparison(original_9km, downsampled_9km, errors, var_name):
    """
    绘制简化比较图表
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 10))

    # 原始数据
    im1 = axes[0].imshow(original_9km, cmap='viridis', origin='lower')
    axes[0].set_title(f'{var_name} - orig-9km')
    #plt.colorbar(im1, ax=axes[0])

    # 降采样数据
    im2 = axes[1].imshow(downsampled_9km, cmap='viridis', origin='lower')
    axes[1].set_title(f'{var_name} - 1km-to-9km')
    #plt.colorbar(im2, ax=axes[1])

    # 误差分布
    im3 = axes[2].imshow(errors, cmap='RdBu_r', origin='lower')
    axes[2].set_title(f'{var_name} - difference')
    #plt.colorbar(im3, ax=axes[2])
    plt.tight_layout()
    plt.show()


def calculate_basic_metrics(era5_25km, downscaled_25km):
    """
    计算基本评估指标：RMSE、MAE、相关系数

    参数:
    - era5_25km: 原始ERA5 25km数据
    - downscaled_25km: 降尺度后的25km数据

    返回:
    - rmse: 均方根误差
    - mae: 平均绝对误差
    - correlation: 相关系数
    - n_points: 有效数据点数
    """

    # 展平数据并移除NaN值
    era5_flat = era5_25km.ravel()
    downscaled_flat = downscaled_25km.ravel()

    n_points = len(era5_flat)

    # 计算误差
    errors = downscaled_flat - era5_flat

    # 计算评估指标
    rmse = np.sqrt(np.mean(errors ** 2))
    mae = np.mean(np.abs(errors))
    correlation, _ = pearsonr(era5_flat, downscaled_flat)

    return rmse, mae, correlation


def print_basic_metrics(rmse, mae, correlation, var_name):
    """
    打印基本评估指标
    """
    print(f"\n=== {var_name} 评估结果 ===")
    print(f"RMSE: {rmse:.4f}")
    print(f"MAE: {mae:.4f}")
    print(f"相关系数: {correlation:.4f}")


def plot_three_images(data1, data2, data3,data4,titles=None, cmap='RdBu_r', figsize=(18, 6)):
    """
    绘制三张图像并排显示

    参数:
    - data1, data2, data3: 三个要显示的数据数组
    - titles: 可选的标题列表 [title1, title2, title3]
    - cmap: 颜色映射，默认为 'RdBu_r'
    - figsize: 图形尺寸，默认为 (18, 6)
    """
    # 设置默认标题
    if titles is None:
        titles = ['orig','Method 1', 'Method 2', 'Method 3']

    # 创建子图
    fig, axs = plt.subplots(1, 4, figsize=figsize)

    # 绘制第一张图
    im0 = axs[0].imshow(np.flipud(data1), cmap=cmap)
    axs[0].set_title(titles[0])
    # 绘制第二张图
    im1 = axs[1].imshow(np.flipud(data2), cmap=cmap)
    axs[1].set_title(titles[1])
    #fig.colorbar(im1, ax=axs[0], location='right')
    # 绘制第三张图
    im2 = axs[2].imshow(np.flipud(data3), cmap=cmap)
    axs[2].set_title(titles[2])
    #fig.colorbar(im2, ax=axs[1], location='right')
    # 绘制第四张图
    im3 = axs[3].imshow(np.flipud(data4), cmap=cmap)
    axs[3].set_title(titles[3])
    #fig.colorbar(im3, ax=axs[2], location='right')

    # 调整布局
    plt.tight_layout()
    plt.show()

def main():
    # 处理数据文件
    compressed_file = "2025-06-01T00_00_00_cn_flatted.nc"
    extracted_data,lat,lon=tiqvdata(load_netcdf_data(compressed_file))
    era5data=nc.Dataset('data_stream-oper_stepType-instant.nc')
    lat_era5 = era5data['latitude'][24:85]
    lon_era5 = era5data['longitude'][4:97]
    #'ERA5',lat_era5[24:85],lon_era5[4:97])
    #由于分辨率导致内存应用不足，将经纬度手动裁切到新疆地区
    lon=lon[11:268]
    lat=lat[189:357]
    #lon[11:268],lat[189:357]
    vars=['t2m', 'v10', 'u10', 'tp6h', 'ssr6h']
    for var in vars:
     analyse_extract_meteorological_data(extracted_data,lat,lon,var=var,era5=era5data[var][:,24:85,4:97],elon=lon_era5,elat=lat_era5)


if __name__ == "__main__":
    main()