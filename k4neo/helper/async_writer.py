import threading
import queue
from k4neo.helper.helper import DiskIO
import pathlib
from loguru import logger


class AsyncDFWriter(threading.Thread):
    def __init__(self, output_path: pathlib.Path, compression: bool = False, max_queue_size=10):
        super().__init__()
        self.queue = queue.Queue(maxsize=max_queue_size)
        self._stop_signal = object()
        self.output_path = output_path
        self.compression = compression
        self.daemon = True

    def run(self):
        while True:
            item = self.queue.get()
            if item is self._stop_signal:
                self.queue.task_done()
                break
            try:
                df, columns, append, header = item
                # logger.info(f"Writing {len(df)} rows to {self.output_path}")
                DiskIO.write_df(
                    df[columns], self.output_path, self.compression, append=append, header=header
                )
            except Exception as e:
                logger.error(f"Failed writing dataframe chunk: {e}")
            finally:
                self.queue.task_done()

    def write(self, df, columns, append, header):
        self.queue.put((df, columns, append, header))

    def wait_until_done(self):
        self.queue.join()

    def stop(self):
        self.queue.put(self._stop_signal)
        self.join()
