import requests
import os
import click
import zipfile


class ZenodoUtils(object):

    def __init__(self, access_token):
        """
        Constructor for ZenodoUtils.
        """
        self._ACCESS_TOKEN = access_token

    def Make_deposition(self, id, path, dir_name):
        url = "https://zenodo.org/api/deposit/depositions/%s?access_token=%s" % (id, self._ACCESS_TOKEN)
        headers = {"Content-Type": "application/json"}

        r = requests.put(url, headers=headers)

        r = requests.post('https://zenodo.org/api/deposit/depositions/%s/actions/publish' % id,
                          params={'access_token': self._ACCESS_TOKEN})

        deposition_id = id
        data = {'name': dir_name}
        files = {'file': open(path, 'rb')}
        r = requests.post('https://zenodo.org/api/deposit/depositions/%s/files' % deposition_id,
                          params={'access_token': self._ACCESS_TOKEN}, data=data,
                          files=files)
        r.json()


def zipdir(path, ziph):
    for root, dirs, files in os.walk(path):
        for file in files:
            ziph.write(os.path.join(root, file),
                       os.path.relpath(os.path.join(root, file),
                                       os.path.join(path, '..')))

@click.command()
@click.option("--id", "-i", type=str, help="Please enter the id")
@click.option("--token", "-t", type=str, help="Please enter the token")
def script(id,token):
    # id = "5189336"
    # token = "DvYIANSfbiDj7CKptZT706vMw2zHdvI1v7Ldn0vTZIY8wSANdDrZVyea3nwb"

    zen = ZenodoUtils(token)

    directory = os.path.join('.', 'SpeciesBed') #enter the path of the directory that holds bed for all the species

    for dir in os.listdir(directory):  # For each folder in the chosen directory
        if dir[0] == '.':  # skipping any hidden file if present
            continue

        path = os.path.join(directory, '%s.zip' % dir)  # path of the file

        zipf = zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED)  # convert the file into zip file by passing the path
        zipdir(dir, zipf)
        zipf.close()
        zen.Make_deposition(id, path, dir)


if __name__ == '__main__':
    script()
