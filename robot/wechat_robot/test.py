from chatterbot import ChatBot
from chatterbot.trainers import ChatterBotCorpusTrainer
import os

os.environ["http_proxy"] = "http://127.0.0.1:7890"
os.environ["https_proxy"] = "http://127.0.0.1:7890"

bot = ChatBot('Sakura',
              storage_adapter='chatterbot.storage.MongoDatabaseAdapter')
trainer = ChatterBotCorpusTrainer(bot)
trainer.train("chatterbot.corpus.chinese")
trainer.train("chatterbot.corpus.english")