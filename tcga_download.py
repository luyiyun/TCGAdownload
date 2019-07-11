# coding:utf-8
'''
这个工具用于下载tcga数据，其使用manifest文件进行下载。
其有断点续传的功能，相对于gdc官方下载工具，其主要的区别是可以在由于网络关系发生
错误从而中断下载后重新发出请求并在断点基础上继续下载，这在下载TCGA大数据（比如
病理图片）时非常有帮助。
此工具配合ss/ssr食用口味更佳。
参数：
    -m manifest文件
    -s 保存的路径
'''
import os
import requests
import argparse
import signal
import time
import datetime
from collections import Counter

import numpy as np
import pandas as pd
import progressbar as pb

print(__doc__)

# 屏蔽warnings
requests.packages.urllib3.disable_warnings()


def get_progressbar(prefix, total_size, start_size=0):
    '''
    得到需要的进度条对象。
    args：
        prefix：进度条前面的信息，string；
        total_size：文件总大小；
        start_size：文件已经下载的大小，用于断点续传；
    return：
        ProgressBar对象
    '''
    widgets = [
        pb.Percentage(), " ",
        pb.Bar(">"), " ",
        pb.AdaptiveETA(), " ",
        pb.AdaptiveTransferSpeed(samples=datetime.timedelta(seconds=1))
    ]
    if prefix is not None:
        widgets = [prefix, " "] + widgets
    bar = pb.ProgressBar(
        widgets=widgets, min_value=start_size, max_value=total_size)
    return bar


class TCGADownload:
    def __init__(self, url, file_path, prefix, total_size=np.inf, sleep_time=5):
        self.url = url
        self.file_path = file_path
        self.prefix = prefix
        self.start_size = 0  # 进度条开始的位置
        self.downloaded_size = 0  # 记录已经下载的量，用于更新进度条
        self.total_size = total_size  # 文件的总大小
        self.sleep_time = sleep_time

    def download(self):
        while self.total_size > self.downloaded_size:
            # 可能网络问题导致提前连接结束，这时需要再次进行请求并断点续传
            download_indx = self._request_part()
            if download_indx:
                self._download_part()
        print('%s have been downloaded' % self.file_path)

    def _request_part(self):
        '''
        return：
            返回一个bool向量，表示是否还要进行download。
        '''
        if self.total_size == np.inf:
            # 如果文件大小未知，则必须进行一次请求去查看文件大小
            self._request()
            self.total_size = int(self.r.headers['content-length'])
            self.open_type = 'wb'  # 如果文件不存在则就会是wb
        if os.path.exists(self.file_path):
            # 如果此文件已经存在，则查看其文件大小
            # ??下面在进行下载的时候，会在最后查看一下downloaded_size和total_size
            #   的大小，如果其不一样大会重新进行请求并断点续传，但这里会再读取
            #   一次文件的大小，而不会充分利用downloaded。。。。
            self.downloaded_size = self.start_size = \
                os.path.getsize(self.file_path)
            if self.downloaded_size == self.total_size:
                return False
            elif self.downloaded_size > self.total_size:
                raise ValueError('现有文件大小大于要下载的文件')
            else:
                headers = {"Range": 'bytes=%d-' % self.start_size}
                self._request(headers=headers)
                self.open_type = 'ab'
        elif self.total_size < np.inf:
            self._request()
            total_size1 = int(self.r.headers['content-length'])
            if total_size1 != self.total_size:
                raise ValueError(
                    '当前已知的文件大小是%d，但网上得到的文件大小是%d。' 
                    % (self.total_size, total_size1)
                )
            self.open_type = 'wb'  # 如果文件不存在则就会是wb
        return True

    def _download_part(self):
        self._progressbar()
        self.bar.start()
        with open(self.file_path, self.open_type) as f:
            for chunk in self.r.iter_content(chunk_size=1024):
                if chunk:
                    self.downloaded_size += len(chunk)
                    f.write(chunk)
                    f.flush()
                self.bar.update(self.downloaded_size)
        self.bar.finish()
        self.r.close()

    def _request(self, headers=None):
        retry = True
        try_errors = []
        while retry:
            try:
                self._close_request()  # 每次在试图请求前都保证之前的请求是关闭的
                self.r = requests.get(
                    self.url, stream=True, verify=False, headers=headers)
                retry = False
            except Exception as e:
                try_errors.append(e.__class__.__name__)
                print(
                    'error:%d, 连接被拒绝，犯的错误是：%s......让我休息5秒钟......'
                    % (len(try_errors), try_errors[-1]), end='\r'
                )
                time.sleep(self.sleep_time)
                continue
        print()  # 前面使用的是\r，这里使得之后的打印可以换到下一行中
        self.request_try_errors = Counter(try_errors)  # 将记录的所有错误都整理成dict
        return self.r

    def _progressbar(self):
        prefix_error = ','.join([
            "%s:%d" % (k, v) for k, v in self.request_try_errors.items()])
        prefixs = self.prefix + ' ' + prefix_error
        self.bar = get_progressbar(prefixs, self.total_size, start_size=0)
        return self.bar

    def _close_request(self):
        if hasattr(self, 'r'):
            self.r.close()


def manifest2df(manifest_path, link):
    df = pd.read_csv(manifest_path, sep='\t', encoding='utf-8')
    if not link.endswith('/'):
        link += '/'
    df['url'] = [link + i for i in df['id'].values]
    return df


def quit(signum, frame):
    # Ctrl+C quit
    print('You choose to stop me.')
    exit()
    print()


def main():
    # Ctrl+C quit
    signal.signal(signal.SIGINT, quit)
    signal.signal(signal.SIGTERM, quit)
    # command line
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m", "--manifest", dest="M", type=str, default="gdc_manifest.txt",
        help="gdc_manifest.txt file path"
    )
    parser.add_argument(
        "-s", "--save", dest="S", type=str, default=os.curdir,
        help="Which folder is the download file saved to?"
    )
    parser.add_argument(
        "-st", "--sleep_time", type=int, default=5,
        help="The time of waiting to next when error occur"
    )
    args = parser.parse_args()

    # args
    link = r'https://api.gdc.cancer.gov/data/'
    manifest_path = args.M
    save_path = args.S
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    print("Save file to {}".format(save_path))

    uuid_df = manifest2df(manifest_path, link=link)
    for i, index in enumerate(uuid_df.index):
        ser = uuid_df.loc[index, :]
        url = ser['url']
        total_size = ser['size']
        file_path = os.path.join(save_path, ser['filename'])
        one_download = TCGADownload(
            url, file_path, total_size=total_size,
            prefix="index: %d, filename: %s" % (i, ser['filename'][:-4]),
            sleep_time=args.sleep_time
        )
        one_download.download()


if __name__ == '__main__':
    main()
