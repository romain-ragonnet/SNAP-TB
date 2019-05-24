import model
import importData

data = importData.data('Indonesia', None)

m = model.TbModel(data, 'scenario_1', -1, initialised=False)
m.store_me()