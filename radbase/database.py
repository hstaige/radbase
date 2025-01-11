import pandas as pd

from radbase.radius_analyzer import RadiiInformation, RadiusDataGateway


class AngeliDataMapper(RadiusDataGateway):

    def __init__(self, file_loc=None):
        self.file_loc = file_loc

    def load_data(self, file_loc=None, bad_types=None) -> RadiiInformation:
        if file_loc is None:
            if self.file_loc is None:
                raise ValueError('file_loc and self.file_loc are both None')
            else:
                file_loc = self.file_loc

        if bad_types is None:
            bad_types = []

        rad_info = RadiiInformation()

        df = pd.read_csv(file_loc)
        for (i, row) in df.iterrows():
            if row['Include'] == 'No' or row['Table'] in bad_types:
                continue
            rad_info.add(term=row['Term'], value=row['Value'], unc=row['Unc'])

            rad_info[max(rad_info.keys())].measurement_type = row['Table']  # TODO: REMOVE THIS HACK

        return rad_info

    def load_analysis(self, **kwargs) -> RadiiInformation:
        pass
