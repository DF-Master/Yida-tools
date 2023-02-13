import telegram

bot = telegram.Bot(token="5927511858:AAHjEsicZsuMJf2pUXsVgOiU0aK-_G43SuM")

for update in bot.getUpdates():
    user_id = update.message.from_user.id
    print(user_id)